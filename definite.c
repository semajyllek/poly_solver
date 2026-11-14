// Definite integrals of rational functions:
//   1. A real root of the denominator in [a, b] means the integral diverges (decided exactly).
//   2. F = integrate(f); the value is F(b) - F(a) with x replaced exactly and constants folded
//      (ln(3/2), atan(1) = π/4, like ln and π terms merged).
//   3. LRT can give atan(N/D) with D vanishing inside (a, b). atan jumps by ±π there, so each
//      crossing adds ∓π times the term's coefficient, making the antiderivative continuous.
//   4. An infinite bound needs deg den >= deg num + 2 and no real pole. Each term's limit comes
//      from its leading behaviour: atan(u) → ±π/2 when u grows, and c·ln|u| → c·ln|lc(u)|
//      (the c·deg(u)·ln|x| parts cancel between the logs, because the integral converges).
//   5. The exact value is compared with adaptive Simpson quadrature of f. On disagreement, or
//      when F has a RootSum (complex logs with branch cuts), only the number is reported.
#include "definite.h"
#include "constfold.h"
#include "roots.h"
#include <math.h>
#include <stdlib.h>

#define AGREE 1e-8 // relative agreement required between exact value and quadrature

// --- real roots in an open interval ---

typedef struct {
    Rational lo, hi;     // lo < hi on the finite sides (0 on an infinite one)
    bool lo_inf, hi_inf; // the span reaches -∞ / +∞
} Span;

// P / gcd(P, P'): the same roots, all simple
static Polynomial* square_free_part(const Polynomial* P) {
    Polynomial* dp = poly_derivative(P);
    Polynomial* g = poly_gcd(P, dp);
    Polynomial* s = poly_divide(P, g);
    free_polynomial(dp);
    free_polynomial(g);
    return s;
}

static bool zero_at(const Polynomial* P, const Rational* x) {
    Rational v;
    rat_init(&v);
    poly_eval_rat(&v, P, x);
    bool z = rat_is_zero(&v);
    rat_clear(&v);
    return z;
}

// P / (x - e) if e is a root of P, else P itself; takes P
static Polynomial* drop_root(Polynomial* P, const Rational* e) {
    if (!zero_at(P, e)) return P;
    Polynomial* lin = create_polynomial();
    Rational c = rat_one_val();
    add_term(lin, &c, 1);
    rat_neg(&c, e);
    add_term(lin, &c, 0);
    Polynomial* q = poly_divide(P, lin);
    rat_clear(&c);
    free_polynomial(lin);
    free_polynomial(P);
    return q;
}

// shrink r until its root (never an endpoint of s) is clearly inside or outside s
static bool root_inside(RootInterval* r, const Polynomial* P, const Span* s) {
    Rational w, two = rat_from_int(2);
    rat_init(&w);
    bool inside;
    for (;;) {
        bool below = !s->lo_inf && rat_cmp(&r->hi, &s->lo) <= 0, above = !s->hi_inf && rat_cmp(&r->lo, &s->hi) >= 0;
        if (below || above) { inside = false; break; }
        bool from = s->lo_inf || rat_cmp(&r->lo, &s->lo) >= 0, to = s->hi_inf || rat_cmp(&r->hi, &s->hi) <= 0;
        if (from && to) { inside = true; break; }
        rat_sub(&w, &r->hi, &r->lo);
        rat_div(&w, &w, &two);
        refine_root(r, P, &w);
    }
    rat_clear(&w);
    rat_clear(&two);
    return inside;
}

// the real roots of P strictly inside s, each isolated in an interval of width <= 1e-9
static RootIntervalList roots_in(const Polynomial* P, const Span* s) {
    RootIntervalList in = { NULL, 0 };
    if (polynomial_degree(P) <= 0) return in;
    Polynomial* sf = square_free_part(P);                 // a root at a finite end is not inside,
    if (!s->lo_inf) sf = drop_root(sf, &s->lo);           // and bisection may never separate it
    if (!s->hi_inf) sf = drop_root(sf, &s->hi);
    RootIntervalList all = isolate_real_roots(sf);
    in.items = malloc((size_t)(all.count + 1) * sizeof(RootInterval));
    Rational w = rat_from_si(1, 1000000000);
    for (int i = 0; i < all.count; i++) {
        if (!root_inside(&all.items[i], sf, s)) continue;
        refine_root(&all.items[i], sf, &w);
        rat_init_set(&in.items[in.count].lo, &all.items[i].lo);
        rat_init_set(&in.items[in.count].hi, &all.items[i].hi);
        in.count++;
    }
    rat_clear(&w);
    root_interval_list_free(&all);
    free_polynomial(sf);
    return in;
}

// P has a root in the closed span
static bool vanishes_on(const Polynomial* P, const Span* s) {
    if ((!s->lo_inf && zero_at(P, &s->lo)) || (!s->hi_inf && zero_at(P, &s->hi))) return true;
    RootIntervalList in = roots_in(P, s);
    bool any = in.count > 0;
    root_interval_list_free(&in);
    return any;
}

// --- numeric value: adaptive Simpson ---

typedef struct {
    const Polynomial* num;
    const Polynomial* den;
} RatFn;

typedef struct {
    double a, b, whole; // Simpson's estimate on [a, b]
} Seg;

typedef struct {
    double eps;   // absolute tolerance for this piece
    int depth;    // remaining bisections
    long* budget; // remaining function evaluations, shared by every piece
} Tol;

// f near an origin o as f(o + u), polynomials in u shifted exactly (so u is small and
// evaluates without cancellation next to a narrow peak at o), on a finite u-range (side 0)
// or on a half-line mapped onto t ∈ [0, 1] by u = side·(1 − t)/t, du = dt/t²
typedef struct {
    RatFn f;     // shifted to the origin
    int side;    // 0, or ±1 for u ∈ [0, +∞) / (-∞, 0]
    double tail; // the mapped integrand at t = 0: lim x²·f(x) as x → ±∞
} Integrand;

static double f_at(const RatFn* f, double x) { return poly_evaluate(f->num, x) / poly_evaluate(f->den, x); }

static double g_at(const Integrand* g, double t) {
    if (!g->side) return f_at(&g->f, t);
    if (t == 0) return g->tail;
    double v = f_at(&g->f, g->side * (1 - t) / t) / (t * t);
    return isfinite(v) ? v : g->tail;                     // overflow far out, where g ≈ its limit
}

static Seg seg_make(const Integrand* f, double a, double b) {
    return (Seg){ a, b, (b - a) / 6 * (g_at(f, a) + 4 * g_at(f, (a + b) / 2) + g_at(f, b)) };
}

// Simpson on both halves, recursing until they agree (with Richardson's correction); the
// depth cap and the evaluation budget bound the work even when round-off stalls convergence
static double asr(const Integrand* f, const Seg* s, Tol t) {
    double m = (s->a + s->b) / 2;
    Seg L = seg_make(f, s->a, m), R = seg_make(f, m, s->b);
    *t.budget -= 6;
    double delta = L.whole + R.whole - s->whole;
    if (t.depth <= 0 || *t.budget <= 0 || fabs(delta) <= 15 * t.eps) return L.whole + R.whole + delta / 15;
    Tol half = { t.eps / 2, t.depth - 1, t.budget };
    return asr(f, &L, half) + asr(f, &R, half);
}

#define PIECES 64 // starting pieces per window, so narrow peaks are seen

// Σ |Simpson| over the starting pieces of [a, b]: the scale the tolerance is relative to
static double rough_of(const Integrand* g, double a, double b) {
    double h = (b - a) / PIECES, rough = 0;
    for (int i = 0; i < PIECES; i++) rough += fabs(seg_make(g, a + i * h, a + (i + 1) * h).whole);
    return rough;
}

// ∫_a^b g, each starting piece refined adaptively to its share of the absolute tolerance t.eps
static double adaptive(const Integrand* g, double a, double b, Tol t) {
    double h = (b - a) / PIECES, total = 0;
    for (int i = 0; i < PIECES; i++) {
        Seg s = seg_make(g, a + i * h, a + (i + 1) * h);
        total += asr(g, &s, (Tol){ t.eps / PIECES, t.depth, t.budget });
    }
    return total;
}

// lim x²·f(x) at ±∞ (the same at both ends): lc(num)/lc(den) when the degrees differ by 2
static double tail_of(const RatFn* f) {
    if (polynomial_degree(f->den) - polynomial_degree(f->num) != 2) return 0;
    Rational n, d;
    rat_init(&n); rat_init(&d);
    poly_coeff(&n, f->num, polynomial_degree(f->num));
    poly_coeff(&d, f->den, polynomial_degree(f->den));
    double t = rat_to_double(&n) / rat_to_double(&d);
    rat_clear(&n); rat_clear(&d);
    return t;
}

// p(x + o) by Horner's rule
static Polynomial* shift(const Polynomial* p, const Rational* o) {
    Polynomial *lin = create_polynomial(), *r = create_polynomial();
    Rational c = rat_one_val();
    add_term(lin, &c, 1);
    add_term(lin, o, 0);
    for (int d = polynomial_degree(p); d >= 0; d--) {
        Polynomial* t = poly_multiply(r, lin);
        poly_coeff(&c, p, d);
        add_term(t, &c, 0);
        free_polynomial(r);
        r = t;
    }
    rat_clear(&c);
    free_polynomial(lin);
    return r;
}

// one quadrature window: f shifted to an origin, over a range of u (t ∈ [0, 1] when mapped)
typedef struct {
    Polynomial *num, *den;
    Integrand g;
    double a, b;
} Window;

// the windows covering the span, for f
typedef struct {
    const RatFn* f;
    Window* w;
    int n;
} Windows;

// u ∈ [u0, u1] (side 0), or the half-line u ≥ 0 (side 1) / u ≤ 0 (side -1)
typedef struct {
    double u0, u1;
    int side;
} Range;

// add the window at origin o over the range r
static void add_window(Windows* W, const Rational* o, Range r) {
    W->w = realloc(W->w, (size_t)(W->n + 1) * sizeof(Window));
    Window* w = &W->w[W->n++];
    w->num = shift(W->f->num, o);
    w->den = shift(W->f->den, o);
    w->g = (Integrand){ { w->num, w->den }, r.side, tail_of(W->f) };
    w->a = r.side ? 0 : r.u0;                             // a half-line is mapped onto t ∈ [0, 1]
    w->b = r.side ? 1 : r.u1;
}

// the windows from p to q, where a NULL p is -∞ and a NULL q is +∞: a half-line anchored at
// its finite end, the whole line at 0, and a finite piece in two halves, each at its own end
static void add_piece(Windows* W, const Rational* p, const Rational* q) {
    Rational zero = rat_zero_val(), h;
    rat_init(&h);
    if (!p && !q) { add_window(W, &zero, (Range){ 0, 0, -1 }); add_window(W, &zero, (Range){ 0, 0, 1 }); }
    else if (!p) add_window(W, q, (Range){ 0, 0, -1 });
    else if (!q) add_window(W, p, (Range){ 0, 0, 1 });
    else {
        rat_sub(&h, q, p);
        double half = rat_to_double(&h) / 2;
        add_window(W, p, (Range){ 0, half, 0 });
        add_window(W, q, (Range){ -half, 0, 0 });
    }
    rat_clear(&zero); rat_clear(&h);
}

// ∫ f over the span, cut at the real critical points of den inside it: a complex pole near
// the real axis makes a narrow peak there that uniform pieces could step over (or the
// half-line map squeeze to nothing), and past the outermost cut |den| is monotone. The
// tolerance is relative 1e-10 of the whole, and the evaluation budget is shared.
static double numeric_value(const RatFn* f, const Span* s) {
    Polynomial* dp = poly_derivative(f->den);
    RootIntervalList cuts = roots_in(dp, s);                   // ascending, width <= 1e-9
    Windows W = { f, NULL, 0 };
    const Rational* p = s->lo_inf ? NULL : &s->lo;
    for (int i = 0; i <= cuts.count; i++) {
        const Rational* q = i < cuts.count ? &cuts.items[i].lo : (s->hi_inf ? NULL : &s->hi);
        add_piece(&W, p, q);
        p = q;
    }
    double rough = 0, v = 0;
    for (int i = 0; i < W.n; i++) rough += rough_of(&W.w[i].g, W.w[i].a, W.w[i].b);
    long budget = 4000000;
    for (int i = 0; i < W.n; i++) {
        v += adaptive(&W.w[i].g, W.w[i].a, W.w[i].b, (Tol){ 1e-10 * (rough + 1e-300) / W.n, 40, &budget });
        free_polynomial(W.w[i].num); free_polynomial(W.w[i].den);
    }
    free(W.w);
    root_interval_list_free(&cuts);
    free_polynomial(dp);
    return v;
}

// --- antiderivative terms and atan jumps ---

typedef struct {
    Expr** t;
    int n;
} Terms;

static void push(Terms* T, Expr* e) {
    T->t = realloc(T->t, (size_t)(T->n + 1) * sizeof(Expr*));
    T->t[T->n++] = e;
}

static Expr* negated(Expr* e) { return create_mul(create_const(-1), e); }

// the summands of e, copied, with their signs applied
static void collect(const Expr* e, int sign, Terms* T) {
    bool add = e->type == EXPR_BINOP && e->binop.op == OP_ADD, sub = e->type == EXPR_BINOP && e->binop.op == OP_SUB;
    if (add || sub) {
        collect(e->binop.left, sign, T);
        collect(e->binop.right, sub ? -sign : sign, T);
    } else if (e->type == EXPR_NEG) {
        collect(e->operand, -sign, T);
    } else {
        push(T, sign > 0 ? expr_copy(e) : negated(expr_copy(e)));
    }
}

// the factor of a product term with the given type (EXPR_ATAN, EXPR_LN), if any
static const Expr* find_fn(const Expr* e, ExprType type) {
    if (e->type == type) return e;
    if (e->type != EXPR_BINOP || e->binop.op != OP_MUL) return NULL;
    const Expr* l = find_fn(e->binop.left, type);
    return l ? l : find_fn(e->binop.right, type);
}

static const Expr* find_atan(const Expr* e) { return find_fn(e, EXPR_ATAN); }

// the nonconstant x-polynomial in a product, if any
static const Polynomial* find_poly(const Expr* e) {
    if (e->type == EXPR_POLY) return polynomial_degree(e->poly) > 0 ? e->poly : NULL;
    if (e->type != EXPR_BINOP || e->binop.op != OP_MUL) return NULL;
    const Polynomial* l = find_poly(e->binop.left);
    return l ? l : find_poly(e->binop.right);
}

static int sign_of(double v) { return (v > 0) - (v < 0); }

// the polynomial W of an atan argument N/(c·W), if it has one
static const Polynomial* atan_pole_poly(const Expr* arg) {
    if (arg->type != EXPR_BINOP || arg->binop.op != OP_DIV) return NULL;
    return find_poly(arg->binop.right);
}

// sign of arg at x (in double precision: only the sign is used)
static int sign_at(const Expr* arg, double x) {
    bool ok = true;
    double v = expr_eval(arg, x, &ok);
    return ok ? sign_of(v) : 0;
}

#define NUDGE 1e-7 // offset from a pole that is an exact rational, for reading signs

// Σ over the poles of arg inside s of the jump of atan(arg), in units of π
static int atan_jumps(const Expr* arg, const Span* s) {
    const Polynomial* W = atan_pole_poly(arg);
    if (!W) return 0;
    RootIntervalList in = roots_in(W, s);
    int jumps = 0;
    for (int i = 0; i < in.count; i++) {                  // ±∞ on either side: ±π/2 → ∓π/2
        double lo = rat_to_double(&in.items[i].lo), hi = rat_to_double(&in.items[i].hi);
        if (rat_eq(&in.items[i].lo, &in.items[i].hi)) { lo -= NUDGE; hi += NUDGE; }
        jumps += (sign_at(arg, hi) - sign_at(arg, lo)) / 2;
    }
    root_interval_list_free(&in);
    return jumps;
}

// copy of the product e with the factor `target` replaced by `with` (taken)
static Expr* replaced(const Expr* e, const Expr* target, Expr* with) {
    if (e == target) return with;
    bool left = find_fn(e->binop.left, target->type) == target;
    return left ? create_mul(replaced(e->binop.left, target, with), expr_copy(e->binop.right))
                : create_mul(expr_copy(e->binop.left), replaced(e->binop.right, target, with));
}

// e with no x left, constants folded; takes e
static Expr* folded(Expr* e) {
    Rational zero = rat_zero_val();
    Expr* v = expr_at(e, &zero);
    free_expr(e);
    rat_clear(&zero);
    return v;
}

// k·π times the coefficient of the atan in term
static Expr* pi_multiple(const Expr* term, const Expr* atan_node, const Rational* k) {
    Expr* c = folded(replaced(term, atan_node, create_const(1)));
    return create_mul(create_const_rat(k), create_mul(c, create_pi()));
}

// --- limits at ±∞ ---

// lim of c·atan(u) or c·ln|u| (dropping c·deg(u)·ln|x|) as x → side·∞, from u ~ L·x^deg
static Expr* fn_limit(const Expr* term, const Expr* node, int side) {
    bool is_atan = node->type == EXPR_ATAN;
    int deg;
    Expr* L = expr_lead(is_atan ? node->fn_arg : node->ln_arg, &deg);
    if (!L) return NULL;
    if (!is_atan) return folded(replaced(term, node, create_ln(L)));
    if (deg <= 0) {                                       // u → L or 0
        if (deg < 0) { free_expr(L); return create_const(0); }
        return folded(replaced(term, node, create_atan(L)));
    }
    bool ok = true;
    int sign = sign_of(expr_eval(L, 0, &ok)) * (deg % 2 ? side : 1);
    free_expr(L);
    Rational half = rat_from_si(sign, 2);                 // atan(±∞) = ±π/2
    Expr* e = ok && sign ? pi_multiple(term, node, &half) : NULL;
    rat_clear(&half);
    return e;
}

// lim of term as x → side·∞; NULL when its leading behaviour is not found or it grows
static Expr* term_limit(const Expr* term, int side) {
    const Expr* node = find_atan(term);
    if (!node) node = find_fn(term, EXPR_LN);
    if (node) return fn_limit(term, node, side);
    int deg;
    Expr* L = expr_lead(term, &deg);                      // a rational term: → 0 or a constant
    if (L && deg < 0) { free_expr(L); return create_const(0); }
    if (L && deg > 0) { free_expr(L); return NULL; }
    return L;
}

// term at the bound x0; at a pole of its atan argument, the limit ±π/2 from inside the range
typedef struct {
    const Rational* x0;
    int inward; // +1 when the range lies to the right of x0
    int inf;    // ±1: the bound is ±∞ and x0 is unused
} Bound;

static Expr* term_at(const Expr* term, const Bound* bd) {
    if (bd->inf) return term_limit(term, bd->inf);
    const Expr* atan_node = find_atan(term);
    const Polynomial* W = atan_node ? atan_pole_poly(atan_node->fn_arg) : NULL;
    if (!W || !zero_at(W, bd->x0)) return expr_at(term, bd->x0);
    double inside = rat_to_double(bd->x0) + bd->inward * NUDGE;
    Rational half = rat_from_si(sign_at(atan_node->fn_arg, inside), 2);
    Expr* e = pi_multiple(term, atan_node, &half);
    rat_clear(&half);
    return e;
}

// the order of two points of the extended real line
static int end_cmp(const Endpoint* p, const Endpoint* q) {
    if (p->inf != q->inf) return p->inf < q->inf ? -1 : 1;
    return p->inf ? 0 : rat_cmp(p->x, q->x);
}

// F(b) - F(a) plus the atan jump corrections over s, folded; NULL when a limit is not found
static Expr* exact_value(const Expr* F, const Bounds* ab, const Span* s) {
    int dir = end_cmp(&ab->a, &ab->b) < 0 ? 1 : -1;
    Terms T = { NULL, 0 }, parts = { NULL, 0 };
    collect(F, 1, &T);
    Bound at_b = { ab->b.x, -dir, ab->b.inf }, at_a = { ab->a.x, dir, ab->a.inf };
    bool ok = true;
    for (int i = 0; i < T.n; i++) {
        Expr* hi = ok ? term_at(T.t[i], &at_b) : NULL;
        Expr* lo = hi ? term_at(T.t[i], &at_a) : NULL;
        ok = hi && lo;
        if (hi) push(&parts, hi);
        if (lo) push(&parts, negated(lo));
        const Expr* atan_node = find_atan(T.t[i]);
        int k = ok && atan_node ? atan_jumps(atan_node->fn_arg, s) : 0;
        Rational kk = rat_from_int(-(long)dir * k);
        if (k) push(&parts, pi_multiple(T.t[i], atan_node, &kk));
        rat_clear(&kk);
        free_expr(T.t[i]);
    }
    free(T.t);
    if (!ok) for (int i = 0; i < parts.n; i++) free_expr(parts.t[i]);
    Expr* v = ok ? const_fold_sum(parts.t, parts.n) : NULL;
    free(parts.t);
    return v;
}

// --- entry point ---

// one definite integral in progress: the integrand, its bounds, and the range they span
typedef struct {
    Expr* f;
    const Bounds* ab;
    const Span* s;
} Problem;

// the exact value when it agrees with the quadrature, else NULL with r->note set
static void exact_or_note(DefiniteResult* r, const Problem* pb) {
    Expr* F = integrate(pb->f, "x");
    if (!F) { r->note = integrate_error(); return; }
    Expr* v = expr_has_rootsum(F) ? NULL : exact_value(F, pb->ab, pb->s);
    bool ok = true;
    double d = v ? expr_eval(v, 0, &ok) : 0;
    if (expr_has_rootsum(F)) r->note = "the antiderivative has RootSum terms, so only the numeric value is given";
    else if (!v) r->note = "the limit of the antiderivative at infinity was not found exactly";
    else if (ok && fabs(d - r->approx) <= AGREE * (1 + fabs(r->approx))) r->exact = v;
    else r->note = "internal check failed: exact value disagrees with quadrature";
    if (v && !r->exact) free_expr(v);
    free_expr(F);
}

// the span from the lower to the upper bound
static Span span_of(const Bounds* ab) {
    bool up = end_cmp(&ab->a, &ab->b) <= 0;
    const Endpoint *lo = up ? &ab->a : &ab->b, *hi = up ? &ab->b : &ab->a;
    Span s = { .lo_inf = lo->inf != 0, .hi_inf = hi->inf != 0 };
    rat_init(&s.lo); rat_init(&s.hi);
    if (!lo->inf) rat_set(&s.lo, lo->x);
    if (!hi->inf) rat_set(&s.hi, hi->x);
    return s;
}

// on an infinite span: f decays like 1/x² or faster (or is 0), so the integral can converge
static bool decays(const RatFn* f) {
    return poly_is_zero(f->num) || polynomial_degree(f->den) - polynomial_degree(f->num) >= 2;
}

DefiniteResult integrate_definite(Expr* f, const Bounds* ab) {
    DefiniteResult r = { NULL, 0, NULL, NULL };
    Expr* c = canonicalize(f);
    if (!c) { r.error = "only polynomials and rational functions can be integrated"; return r; }
    Polynomial* one = create_polynomial();
    Rational u = rat_one_val();
    add_term(one, &u, 0);
    rat_clear(&u);
    RatFn fn = { c->type == EXPR_POLY ? c->poly : c->ratfn.num, c->type == EXPR_POLY ? one : c->ratfn.den };
    Span s = span_of(ab);
    int dir = end_cmp(&ab->a, &ab->b);
    if ((s.lo_inf || s.hi_inf) && !decays(&fn))            // before the empty range: ∫_∞^∞ x² is not 0
        r.error = "the integral diverges at infinity (the denominator's degree must exceed the numerator's by 2 or more)";
    else if (dir == 0 || poly_is_zero(fn.num)) r.exact = create_const(0);
    else if (vanishes_on(fn.den, &s)) r.error = "the integrand has a pole in [a, b], so the integral diverges";
    else {
        r.approx = (dir < 0 ? 1 : -1) * numeric_value(&fn, &s);
        exact_or_note(&r, &(Problem){ f, ab, &s });
    }
    rat_clear(&s.lo); rat_clear(&s.hi);
    free_polynomial(one);
    free_expr(c);
    return r;
}

void definite_free(DefiniteResult* r) {
    if (r->exact) free_expr(r->exact);
    r->exact = NULL;
}
