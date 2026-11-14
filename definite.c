// Definite integrals of rational functions:
//   1. A real root of the denominator in [a, b] means the integral diverges (decided exactly).
//   2. F = integrate(f); the value is F(b) - F(a) with x replaced exactly and constants folded
//      (ln(3/2), atan(1) = π/4, like ln and π terms merged).
//   3. LRT can give atan(N/D) with D vanishing inside (a, b). atan jumps by ±π there, so each
//      crossing adds ∓π times the term's coefficient, making the antiderivative continuous.
//   4. The exact value is compared with adaptive Simpson quadrature of f. On disagreement, or
//      when F has a RootSum (complex logs with branch cuts), only the number is reported.
#include "definite.h"
#include "constfold.h"
#include "roots.h"
#include <math.h>
#include <stdlib.h>

#define AGREE 1e-8 // relative agreement required between exact value and quadrature

// --- real roots in an open interval ---

typedef struct {
    Rational lo, hi; // lo < hi
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

// shrink r until its root (never an endpoint of s) is clearly inside or outside s
static bool root_inside(RootInterval* r, const Polynomial* P, const Span* s) {
    Rational w, two = rat_from_int(2);
    rat_init(&w);
    bool inside;
    for (;;) {
        if (rat_cmp(&r->hi, &s->lo) <= 0 || rat_cmp(&r->lo, &s->hi) >= 0) { inside = false; break; }
        if (rat_cmp(&r->lo, &s->lo) >= 0 && rat_cmp(&r->hi, &s->hi) <= 0) { inside = true; break; }
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
    Polynomial* sf = square_free_part(P);
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

static bool zero_at(const Polynomial* P, const Rational* x) {
    Rational v;
    rat_init(&v);
    poly_eval_rat(&v, P, x);
    bool z = rat_is_zero(&v);
    rat_clear(&v);
    return z;
}

// P has a root in the closed interval [lo, hi]
static bool vanishes_on(const Polynomial* P, const Span* s) {
    if (zero_at(P, &s->lo) || zero_at(P, &s->hi)) return true;
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

static double f_at(const RatFn* f, double x) { return poly_evaluate(f->num, x) / poly_evaluate(f->den, x); }

static Seg seg_make(const RatFn* f, double a, double b) {
    return (Seg){ a, b, (b - a) / 6 * (f_at(f, a) + 4 * f_at(f, (a + b) / 2) + f_at(f, b)) };
}

// Simpson on both halves, recursing until they agree (with Richardson's correction); the
// depth cap and the evaluation budget bound the work even when round-off stalls convergence
static double asr(const RatFn* f, const Seg* s, Tol t) {
    double m = (s->a + s->b) / 2;
    Seg L = seg_make(f, s->a, m), R = seg_make(f, m, s->b);
    *t.budget -= 6;
    double delta = L.whole + R.whole - s->whole;
    if (t.depth <= 0 || *t.budget <= 0 || fabs(delta) <= 15 * t.eps) return L.whole + R.whole + delta / 15;
    Tol half = { t.eps / 2, t.depth - 1, t.budget };
    return asr(f, &L, half) + asr(f, &R, half);
}

// ∫_a^b f numerically: 64 starting pieces (so narrow peaks are seen), relative tolerance 1e-10
static double quadrature(const RatFn* f, double a, double b) {
    const int pieces = 64;
    double h = (b - a) / pieces, rough = 0, total = 0;
    for (int i = 0; i < pieces; i++) rough += fabs(seg_make(f, a + i * h, a + (i + 1) * h).whole);
    long budget = 4000000;
    for (int i = 0; i < pieces; i++) {
        Seg s = seg_make(f, a + i * h, a + (i + 1) * h);
        total += asr(f, &s, (Tol){ 1e-10 * (rough + 1e-300) / pieces, 40, &budget });
    }
    return total;
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

// the atan factor of a product term, if any
static const Expr* find_atan(const Expr* e) {
    if (e->type == EXPR_ATAN) return e;
    if (e->type != EXPR_BINOP || e->binop.op != OP_MUL) return NULL;
    const Expr* l = find_atan(e->binop.left);
    return l ? l : find_atan(e->binop.right);
}

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

// copy of e with the node `target` replaced by 1
static Expr* without(const Expr* e, const Expr* target) {
    if (e == target) return create_const(1);
    if (e->type == EXPR_BINOP && e->binop.op == OP_MUL)
        return create_mul(without(e->binop.left, target), without(e->binop.right, target));
    return expr_copy(e);
}

// k·π times the coefficient of the atan in term
static Expr* pi_multiple(const Expr* term, const Expr* atan_node, const Rational* k) {
    Rational zero = rat_zero_val();
    Expr* coef = without(term, atan_node);
    Expr* c = expr_at(coef, &zero);                       // no x left: just folds the constants
    free_expr(coef);
    rat_clear(&zero);
    return create_mul(create_const_rat(k), create_mul(c, create_pi()));
}

// term at the bound x0; at a pole of its atan argument, the limit ±π/2 from inside the range
typedef struct {
    const Rational* x0;
    int inward; // +1 when the range lies to the right of x0
} Bound;

static Expr* term_at(const Expr* term, const Bound* bd) {
    const Expr* atan_node = find_atan(term);
    const Polynomial* W = atan_node ? atan_pole_poly(atan_node->fn_arg) : NULL;
    if (!W || !zero_at(W, bd->x0)) return expr_at(term, bd->x0);
    double inside = rat_to_double(bd->x0) + bd->inward * NUDGE;
    Rational half = rat_from_si(sign_at(atan_node->fn_arg, inside), 2);
    Expr* e = pi_multiple(term, atan_node, &half);
    rat_clear(&half);
    return e;
}

// F(b) - F(a) plus the atan jump corrections over s, folded
static Expr* exact_value(const Expr* F, const Bounds* ab, const Span* s) {
    int dir = rat_cmp(ab->a, ab->b) < 0 ? 1 : -1;
    Terms T = { NULL, 0 }, parts = { NULL, 0 };
    collect(F, 1, &T);
    Bound at_b = { ab->b, -dir }, at_a = { ab->a, dir };
    for (int i = 0; i < T.n; i++) {
        push(&parts, term_at(T.t[i], &at_b));
        push(&parts, negated(term_at(T.t[i], &at_a)));
        const Expr* atan_node = find_atan(T.t[i]);
        int k = atan_node ? atan_jumps(atan_node->fn_arg, s) : 0;
        Rational kk = rat_from_int(-(long)dir * k);
        if (k) push(&parts, pi_multiple(T.t[i], atan_node, &kk));
        rat_clear(&kk);
        free_expr(T.t[i]);
    }
    free(T.t);
    Expr* v = const_fold_sum(parts.t, parts.n);
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
    Expr* f = pb->f;
    const Bounds* ab = pb->ab;
    const Span* s = pb->s;
    Expr* F = integrate(f, "x");
    if (!F) { r->note = integrate_error(); return; }
    if (expr_has_rootsum(F)) {
        r->note = "the antiderivative has RootSum terms, so only the numeric value is given";
    } else {
        Expr* v = exact_value(F, ab, s);
        bool ok = true;
        double d = expr_eval(v, 0, &ok);
        if (ok && fabs(d - r->approx) <= AGREE * (1 + fabs(r->approx))) r->exact = v;
        else { free_expr(v); r->note = "internal check failed: exact value disagrees with quadrature"; }
    }
    free_expr(F);
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
    Span s;
    bool up = rat_cmp(ab->a, ab->b) <= 0;
    rat_init_set(&s.lo, up ? ab->a : ab->b);
    rat_init_set(&s.hi, up ? ab->b : ab->a);
    if (rat_eq(&s.lo, &s.hi)) r.exact = create_const(0);
    else if (vanishes_on(fn.den, &s)) r.error = "the integrand has a pole in [a, b], so the integral diverges";
    else {
        r.approx = quadrature(&fn, rat_to_double(ab->a), rat_to_double(ab->b));
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
