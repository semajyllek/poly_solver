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
// Exponential sums f = c_0 + Σ c_w·e^w (exp, sin, cos; expsum.h) and bounds q·π follow the
// same steps, with these differences (definite_exp):
//   * a root of a coefficient denominator in [a, b] is reported as unsupported, not as a pole
//     (sin(x)/x is finite at 0), unless f = c_0 is rational;
//   * an infinite end needs Re w → -∞ there for every w ≠ 0 (then those terms of F tend to 0)
//     and c_0 decaying like 1/x²; slower decay of c_0 is divergence, while an oscillating
//     term such as cos(x)/x is unsupported (quadrature cannot confirm a conditional limit);
//   * q·π is replaced by a 50-digit rational inside a band that must hold no root of the
//     relevant polynomials (a rational polynomial never vanishes at q·π, so the band only has
//     to be narrow enough), and F is evaluated at q·π exactly (constfold.h: expr_at_pi);
//   * without an elementary antiderivative, the numeric value is given alone.
#include "definite.h"
#include "constfold.h"
#include "expsum.h"
#include "deftower.h"
#include "tanfield.h"
#include "risch.h"
#include "roots.h"
#include "parser.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define AGREE 1e-8 // relative agreement required between exact value and quadrature

// --- real roots in an open interval ---

typedef struct {
    Rational lo, hi;     // lo < hi on the finite sides (0 on an infinite one); q·π as a 50-digit rational
    bool lo_inf, hi_inf; // the span reaches -∞ / +∞
    const Bounds* ab;    // the bounds it was made from (for the bands around q·π)
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

// --- bounds q·π ---

// the value of a finite end: x, or x·PI_50
static void end_value(const Endpoint* p, Rational* v) {
    rat_set(v, p->x);
    if (!p->pi) return;
    Rational pi;
    rat_init(&pi);
    pi_50(&pi);
    rat_mul(v, v, &pi);
    rat_clear(&pi);
}

// the band [v - δ, v + δ] around a q·π end, δ = (|q| + 1)·1e-49, holding the true q·π
static Span band_of(const Endpoint* p) {
    Span b = { .lo_inf = false, .hi_inf = false, .ab = NULL };
    rat_init(&b.lo); rat_init(&b.hi);
    Rational v, d;
    rat_init(&v); rat_init(&d);
    end_value(p, &v);
    rat_abs(&d, p->x);
    Rational one = rat_one_val(), eps;
    rat_init(&eps);
    rat_add(&d, &d, &one);
    rat_parse_literal(&eps, "1e-49");
    rat_mul(&d, &d, &eps);
    rat_sub(&b.lo, &v, &d);
    rat_add(&b.hi, &v, &d);
    rat_clear(&v); rat_clear(&d); rat_clear(&one); rat_clear(&eps);
    return b;
}

// no q·π end of s has a root of P (or the other, rational, end) in its band, so the
// 50-digit stand-in for q·π sees the same roots on each side
static bool clear_of_bands(const Polynomial* P, const Span* s) {
    const Endpoint* ends[2] = { &s->ab->a, &s->ab->b };
    bool clear = true;
    for (int i = 0; i < 2 && clear; i++) {
        if (!ends[i]->pi || ends[i]->inf) continue;
        Span b = band_of(ends[i]);
        clear = !vanishes_on(P, &b);
        const Endpoint* o = ends[1 - i];
        if (clear && !o->inf && !o->pi) clear = rat_cmp(o->x, &b.lo) < 0 || rat_cmp(o->x, &b.hi) > 0;
        rat_clear(&b.lo); rat_clear(&b.hi);
    }
    return clear;
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

// a polynomial p near an origin o, in doubles: p(o + u) = c(u + off), with c either p shifted
// exactly to o (off = 0) or p itself (off = o), whichever rounds less over the window
typedef struct {
    DPoly c;
    double off;
} Near;

static double near_at(const Near* p, double u) { return dpoly_eval(&p->c, u + p->off); }

// f near an origin o as f(o + u) (Near: shifted, so u is small and evaluates without
// cancellation next to a narrow peak at o, or not, when the shift itself would cancel), on a
// finite u-range (side 0) or on a half-line mapped onto t ∈ [0, 1] by u = side·(1 − t)/t, du = dt/t²
typedef struct {
    Near n, d;     // f's numerator and denominator
    int side;      // 0, or ±1 for u ∈ [0, +∞) / (-∞, 0]
    double tail;   // the mapped integrand at t = 0: lim x²·f(x) as x → ±∞
    const Expr* e; // instead of f: an expression, evaluated at x = o + u
    double o;
} Integrand;

static double f_at(const Integrand* g, double u) { return near_at(&g->n, u) / near_at(&g->d, u); }

// the integrand at u (from its origin)
static double u_at(const Integrand* g, double u) {
    if (!g->e) return f_at(g, u);
    bool ok = true;
    double v = expr_eval(g->e, g->o + u, &ok);
    return ok ? v : NAN;
}

static double g_at(const Integrand* g, double t) {
    if (!g->side) return u_at(g, t);
    if (t == 0) return g->tail;
    double v = u_at(g, g->side * (1 - t) / t) / (t * t);
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

// one quadrature window: f near an origin, over a range of u (t ∈ [0, 1] when mapped)
typedef struct {
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

#define SHIFT_MAX_DEGREE 32 // the exact shift costs O(deg²) big rationals; above this p stays as it is

// Horner's error bound for c at t, up to the unit round-off: Σ |c_i|·|t|^i
static double horner_bound(const DPoly* c, double t) {
    double b = 0;
    for (int i = c->deg; i >= 0; i--) b = b * fabs(t) + fabs(c->c[i]);
    return b;
}

// p near o over r: shifted when that bounds the round-off lower at the window's far end |u|
// (or u = 1 on a half-line, where the peak next to o matters), else p's own coefficients;
// (1 + u)^100 from x^100 at o = 1 cancels catastrophically, x² - 100x + 2500.0001 at 50 does not
static Near near_of(const Polynomial* p, const Rational* o, const Range* r) {
    double od = rat_to_double(o), u = r->side ? 1 : fmax(fabs(r->u0), fabs(r->u1));
    Near plain = { dpoly_from(p), od };
    if (polynomial_degree(p) > SHIFT_MAX_DEGREE) return plain;
    Polynomial* s = shift(p, o);
    Near shifted = { dpoly_from(s), 0 };
    free_polynomial(s);
    bool use = horner_bound(&shifted.c, u) <= horner_bound(&plain.c, od + u * (r->side ? r->side : 1));
    dpoly_free(use ? &plain.c : &shifted.c);
    return use ? shifted : plain;
}

#define WIDE_WINDOW 64.0  // a finite window wider than this many peak scales is cut geometrically
#define MAX_PIECES 1100   // geometric pieces per window: 2^1100 spans every double

// the window at origin o over the range r
static void push_window(Windows* W, const Rational* o, Range r) {
    W->w = realloc(W->w, (size_t)(W->n + 1) * sizeof(Window));
    Window* w = &W->w[W->n++];
    w->g = (Integrand){ .n = near_of(W->f->num, o, &r), .d = near_of(W->f->den, o, &r), .side = r.side, .tail = tail_of(W->f) };
    w->a = r.side ? 0 : r.u0;                             // a half-line is mapped onto t ∈ [0, 1]
    w->b = r.side ? 1 : r.u1;
}

// the width of a peak of 1/den at o: half the least |c_0/c_k|^(1/k) of den(o + u) = Σ c_k·u^k,
// a lower bound on the distance from o to den's nearest root; 1 above SHIFT_MAX_DEGREE
static double peak_scale(const Polynomial* den, const Rational* o) {
    if (polynomial_degree(den) > SHIFT_MAX_DEGREE) return 1;
    Polynomial* s = shift(den, o);
    DPoly c = dpoly_from(s);
    double m = INFINITY;
    for (int k = 1; k <= c.deg; k++)
        if (c.c[k] != 0) m = fmin(m, pow(fabs(c.c[0] / c.c[k]), 1.0 / k));
    dpoly_free(&c);
    free_polynomial(s);
    return m / 2;
}

// add the window at origin o over r; a finite one, [0, h] or [-h, 0], wider than WIDE_WINDOW peak
// scales s goes in pieces [0, s], [s, 2s], [2s, 4s], ... toward h, so a peak at o is not lost
// among uniform pieces of width h/64 (1/(x^2 + 1) on [0, 1e40], 1/(x^2 + 1e-300) on [0, 1])
static void add_window(Windows* W, const Rational* o, Range r) {
    double h = r.u0 + r.u1;                               // the far end, signed: the other end is 0
    double s = r.side ? INFINITY : fmax(peak_scale(W->f->den, o), fabs(h) * ldexp(1, -MAX_PIECES));
    if (r.side || !(fabs(h) > WIDE_WINDOW * s)) { push_window(W, o, r); return; }
    for (double a = 0, b = s; a < fabs(h); a = b, b *= 2) {
        double e = fmin(b, fabs(h));
        push_window(W, o, h > 0 ? (Range){ a, e, 0 } : (Range){ -e, -a, 0 });
        b = e;
    }
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
        dpoly_free(&W.w[i].g.n.c); dpoly_free(&W.w[i].g.d.c);
    }
    free(W.w);
    root_interval_list_free(&cuts);
    free_polynomial(dp);
    return v;
}

// ∫ e over the span, for an e that is finite on it and decays exponentially at an infinite
// end (tail 0): [lo, hi], or half-lines from the finite end (from 0 for the whole line).
// *scale is the rough ∫|e| the tolerance is relative to
static double numeric_value_expr(const Expr* e, const Span* s, double* scale) {
    double lo = rat_to_double(&s->lo), hi = rat_to_double(&s->hi);
    Integrand g[2];
    double b[2] = { 1, 1 };
    int n = 0;
    if (!s->lo_inf && !s->hi_inf) { g[n] = (Integrand){ .e = e, .o = lo }; b[n++] = hi - lo; }
    else if (!s->lo_inf) g[n++] = (Integrand){ .e = e, .o = lo, .side = 1 };
    else if (!s->hi_inf) g[n++] = (Integrand){ .e = e, .o = hi, .side = -1 };
    else { g[n++] = (Integrand){ .e = e, .side = -1 }; g[n++] = (Integrand){ .e = e, .side = 1 }; }
    double rough = 0, v = 0;
    for (int i = 0; i < n; i++) rough += rough_of(&g[i], 0, b[i]);
    long budget = 4000000;
    for (int i = 0; i < n; i++) v += adaptive(&g[i], 0, b[i], (Tol){ 1e-10 * (rough + 1e-300) / n, 40, &budget });
    *scale = rough;
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
static int atan_jumps(const Expr* arg, const Span* s, bool* ok) {
    const Polynomial* W = atan_pole_poly(arg);
    if (!W) return 0;
    if (s->ab && !clear_of_bands(W, s)) { *ok = false; return 0; }
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

// lim of term as x → side·∞; NULL when its leading behaviour is not found or it grows. A term
// with exp, sin or cos is 0 there: the caller has checked that every e^w decays at that end
static Expr* term_limit(const Expr* term, int side) {
    if (expr_has_exp(term) || expr_has_trig(term)) return create_const(0);
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
    int inward;            // +1 when the range lies to the right of x0
    int inf;               // ±1: the bound is ±∞ and x0 is unused
    bool pi;               // the bound is x0·π (never a pole of an atan argument: π is transcendental)
    const Polynomial* H;   // θ = tan(H) in the answer, or NULL
    bool tan_inf;          // tan(H) is infinite at the bound
} Bound;

static double bound_double(const Bound* bd) { return rat_to_double(bd->x0) * (bd->pi ? M_PI : 1); }

// a polynomial argument equal to H
static bool is_arg(const Expr* a, const Polynomial* H) {
    if (a->type != EXPR_POLY) return false;
    Polynomial* d = poly_subtract(a->poly, H);
    bool same = poly_is_zero(d);
    free_polynomial(d);
    return same;
}

// e contains tan(H), the tree sin(H)/cos(H)
static bool has_tan(const Expr* e, const Polynomial* H) {
    if (!H) return false;
    switch (e->type) {
        case EXPR_BINOP: {
            const Expr *l = e->binop.left, *r = e->binop.right;
            if (e->binop.op == OP_DIV && l->type == EXPR_SIN && r->type == EXPR_COS && is_arg(l->fn_arg, H)) return true;
            return has_tan(l, H) || has_tan(r, H);
        }
        case EXPR_NEG: return has_tan(e->operand, H);
        case EXPR_LN: return has_tan(e->ln_arg, H);
        case EXPR_ATAN: case EXPR_SQRT: case EXPR_EXP: case EXPR_SIN: case EXPR_COS: return has_tan(e->fn_arg, H);
        default: return false;
    }
}

// the term at a bound where tan(H) is infinite: atan(P(tan H)) → ±π/2 from inside, a proper
// fraction in tan(H) → 0, a term without tan(H) as it is; NULL (only the number is then
// given) for a log of tan(H) or of cos(H), whose infinite parts cancel between terms
static Expr* term_at_tan_infinity(const Expr* term, const Bound* bd) {
    const Expr* a = find_atan(term);
    if (a && has_tan(a->fn_arg, bd->H)) {
        Rational half = rat_from_si(sign_at(a->fn_arg, bound_double(bd) + bd->inward * NUDGE), 2);
        Expr* e = pi_multiple(term, a, &half);
        rat_clear(&half);
        return e;
    }
    const Expr* ln = find_fn(term, EXPR_LN);
    if (ln && (has_tan(ln->ln_arg, bd->H) || (ln->ln_arg->type == EXPR_COS && is_arg(ln->ln_arg->fn_arg, bd->H)))) return NULL;
    if (!has_tan(term, bd->H)) return bd->pi ? expr_at_pi(term, bd->x0) : expr_at(term, bd->x0);
    KFrac q;
    if (!tan_frac(term, bd->H, &q)) return NULL;
    bool proper = kp_deg(&q.n) < kp_deg(&q.d);
    kfrac_free(&q);
    return proper ? create_const(0) : NULL;
}

static Expr* term_at(const Expr* term, const Bound* bd) {
    if (bd->inf) return term_limit(term, bd->inf);
    if (bd->tan_inf) return term_at_tan_infinity(term, bd);
    if (bd->pi) return expr_at_pi(term, bd->x0);
    const Expr* atan_node = find_atan(term);
    const Polynomial* W = atan_node ? atan_pole_poly(atan_node->fn_arg) : NULL;
    if (!W || !zero_at(W, bd->x0)) return expr_at(term, bd->x0);
    double inside = rat_to_double(bd->x0) + bd->inward * NUDGE;
    Rational half = rat_from_si(sign_at(atan_node->fn_arg, inside), 2);
    Expr* e = pi_multiple(term, atan_node, &half);
    rat_clear(&half);
    return e;
}

// the order of two points of the extended real line (q·π against a rational by its 50-digit
// value: clear_of_bands rules out a tie inside the band)
static int end_cmp(const Endpoint* p, const Endpoint* q) {
    if (p->inf != q->inf) return p->inf < q->inf ? -1 : 1;
    if (p->inf) return 0;
    if (p->pi == q->pi) return rat_cmp(p->x, q->x);
    Rational u, v;
    rat_init(&u); rat_init(&v);
    end_value(p, &u); end_value(q, &v);
    int c = rat_cmp(&u, &v);
    rat_clear(&u); rat_clear(&v);
    return c;
}

// one definite integral in progress: the integrand, its bounds, the range they span, and for
// an answer in θ = tan(H) the points inside where θ is infinite
typedef struct {
    Expr* f;
    const Bounds* ab;
    const Span* s;
    double scale;            // ∫|f| roughly: agreement is also relative to it (F(b) - F(a) may cancel to 0)
    const Polynomial* H;     // or NULL
    const TowerPoints* tp;   // with H
} Problem;

// the jumps of atan(P(tan H)) where tan(H) is infinite inside the range, in units of π
static int tan_jumps(const Expr* arg, const Problem* pb) {
    if (!pb->tp || !has_tan(arg, pb->H)) return 0;
    int k = 0;
    for (int i = 0; i < pb->tp->n; i++)
        k += (sign_at(arg, pb->tp->cross[i] + NUDGE) - sign_at(arg, pb->tp->cross[i] - NUDGE)) / 2;
    return k;
}

// F(b) - F(a) plus the atan jump corrections over the range, folded; NULL when a limit is not found
static Expr* exact_value(const Expr* F, const Problem* pb) {
    const Bounds* ab = pb->ab;
    int dir = end_cmp(&ab->a, &ab->b) < 0 ? 1 : -1;
    Terms T = { NULL, 0 }, parts = { NULL, 0 };
    collect(F, 1, &T);
    bool b_inf = pb->tp && (dir > 0 ? pb->tp->at_hi : pb->tp->at_lo), a_inf = pb->tp && (dir > 0 ? pb->tp->at_lo : pb->tp->at_hi);
    Bound at_b = { ab->b.x, -dir, ab->b.inf, ab->b.pi, pb->H, b_inf }, at_a = { ab->a.x, dir, ab->a.inf, ab->a.pi, pb->H, a_inf };
    bool ok = true;
    for (int i = 0; i < T.n; i++) {
        Expr* hi = ok ? term_at(T.t[i], &at_b) : NULL;
        Expr* lo = hi ? term_at(T.t[i], &at_a) : NULL;
        ok = hi && lo;
        if (hi) push(&parts, hi);
        if (lo) push(&parts, negated(lo));
        const Expr* atan_node = find_atan(T.t[i]);
        int k = ok && atan_node ? atan_jumps(atan_node->fn_arg, pb->s, &ok) + tan_jumps(atan_node->fn_arg, pb) : 0;
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

// the exact value from the antiderivative F (taken) when it agrees with the quadrature, else
// NULL with r->note set
static void exact_or_note(DefiniteResult* r, const Problem* pb, Expr* F) {
    Expr* v = expr_has_rootsum(F) ? NULL : exact_value(F, pb);
    bool ok = true;
    double d = v ? expr_eval(v, 0, &ok) : 0;
    if (expr_has_rootsum(F)) r->note = "the antiderivative has RootSum terms, so only the numeric value is given";
    else if (!v) r->note = "a limit of the antiderivative (at infinity, or where tan is infinite) was not found exactly";
    else if (ok && fabs(d - r->approx) <= AGREE * (1 + fabs(r->approx) + pb->scale)) r->exact = v;
    else if (ok && isinf(d)) {                            // beyond doubles (x^2 on [0, 1e200]): nothing to compare
        r->exact = v;
        r->approx = d;
        r->note = "beyond double precision, so not compared with quadrature";
    }
    else r->note = "internal check failed: exact value disagrees with quadrature";
    if (v && !r->exact) free_expr(v);
    free_expr(F);
}

// the rational integrator's antiderivative, or NULL with r->note set
static Expr* rational_antiderivative(DefiniteResult* r, const Expr* f) {
    Expr* F = integrate((Expr*)f, "x");
    if (!F) r->note = integrate_error();
    return F;
}

// the span from the lower to the upper bound
static Span span_of(const Bounds* ab) {
    bool up = end_cmp(&ab->a, &ab->b) <= 0;
    const Endpoint *lo = up ? &ab->a : &ab->b, *hi = up ? &ab->b : &ab->a;
    Span s = { .lo_inf = lo->inf != 0, .hi_inf = hi->inf != 0, .ab = ab };
    rat_init(&s.lo); rat_init(&s.hi);
    if (!lo->inf) end_value(lo, &s.lo);
    if (!hi->inf) end_value(hi, &s.hi);
    return s;
}

// on an infinite span: f decays like 1/x² or faster (or is 0), so the integral can converge
static bool decays(const RatFn* f) {
    return poly_is_zero(f->num) || polynomial_degree(f->den) - polynomial_degree(f->num) >= 2;
}

#define DIVERGES_AT_INF "the integral diverges at infinity (the denominator's degree must exceed the numerator's by 2 or more)"

// why ∫ over s cannot be taken in the exponential-sum setting, or NULL
static const char* exp_obstacle(const ExpSumShape* sh, const Span* s) {
    RatFn c0 = { sh->c0_num, sh->c0_den };
    if (!clear_of_bands(sh->poles, s)) return "a bound q·pi lies too close to a pole or to the other bound to be decided";
    const bool inf[2] = { s->lo_inf, s->hi_inf };
    for (int side = 0; side < 2; side++) {
        if (!inf[side]) continue;
        if (!sh->decays[side])
            return "at infinity the integrand must decay exponentially (exp(E) with E -> -inf in every term); oscillating terms such as cos(x)/x are not supported yet";
        if (!decays(&c0)) return DIVERGES_AT_INF;
    }
    if (vanishes_on(sh->poles, s))
        return sh->pure ? "the integrand has a pole in [a, b], so the integral diverges"
                        : "a coefficient's denominator vanishes in [a, b] (as for sin(x)/x at 0); this is not supported yet";
    return NULL;
}

// the antiderivative of an exponential sum, or NULL with r->note set (no elementary one,
// or outside the integrator's reach: then only the number is given)
static Expr* exp_antiderivative(DefiniteResult* r, const Expr* f, const ExpSumShape* sh) {
    if (sh->pure) {                                       // c_0 itself: f may still spell 0·cos(x)
        Expr* c0 = create_rational_fn(poly_copy(sh->c0_num), poly_copy(sh->c0_den));
        Expr* F = rational_antiderivative(r, c0);
        free_expr(c0);
        return F;
    }
    RischResult I = integrate_exp(f);
    if (I.status == RISCH_OK) return I.result;
    free_expr(I.result);
    r->note = I.status == RISCH_NONELEMENTARY ? "no elementary antiderivative, so only the numeric value is given" : I.reason;
    return NULL;
}

// the product of the x-denominators of N's coefficients (f = N/D in θ)
static Polynomial* coefficient_poles(const LPoly* N) {
    Rational one = rat_one_val();
    Polynomial* p = poly_constant(&one);
    rat_clear(&one);
    for (int i = 0; i < N->n; i++) { Polynomial* q = poly_multiply(p, N->a[i].den); free_polynomial(p); p = q; }
    return p;
}

// why the tower form F of f cannot be integrated over s exactly, or NULL (tp filled)
static const char* tower_obstacle(const ExpForm* F, const Span* s, TowerPoints* tp) {
    *tp = (TowerPoints){ NULL, 0, false, false };
    if (F->unsupported) return F->unsupported;
    if (s->lo_inf || s->hi_inf) return "infinite bounds with exp, sin or cos in a denominator are not supported yet";
    Polynomial* P = coefficient_poles(&F->f);
    bool vanishes = vanishes_on(P, s) || !clear_of_bands(P, s);
    free_polynomial(P);
    if (vanishes) return "a coefficient's denominator vanishes in [a, b]; this is not supported yet";
    bool up = end_cmp(&s->ab->a, &s->ab->b) <= 0;
    TowerEnds E = { up ? &s->ab->a : &s->ab->b, up ? &s->ab->b : &s->ab->a };
    return tower_points(F, &E, tp);
}

// f = N/D with a sum of exp or of sin/cos in the denominator, over a finite range
static DefiniteResult definite_tower(Expr* f, const Bounds* ab) {
    DefiniteResult r = { NULL, 0, NULL, NULL };
    ExpForm F = expr_has_trig(f) ? tan_form(f) : exp_form(f);
    Span s = span_of(ab);
    TowerPoints tp;
    r.error = tower_obstacle(&F, &s, &tp);
    if (!r.error && end_cmp(&ab->a, &ab->b) == 0) r.exact = create_const(0);
    else if (!r.error) {
        double scale = 0;
        r.approx = (end_cmp(&ab->a, &ab->b) < 0 ? 1 : -1) * numeric_value_expr(f, &s, &scale);
        RischResult I = integrate_exp(f);
        if (I.status != RISCH_OK) {
            free_expr(I.result);
            r.note = I.status == RISCH_NONELEMENTARY ? "no elementary antiderivative, so only the numeric value is given" : I.reason;
        } else {
            bool tan = F.t.kind == TOWER_TAN;
            exact_or_note(&r, &(Problem){ f, ab, &s, scale, tan ? F.t.v : NULL, tan ? &tp : NULL }, I.result);
        }
    }
    tower_points_free(&tp);
    rat_clear(&s.lo); rat_clear(&s.hi);
    exp_form_free(&F);
    return r;
}

// f with exp, sin or cos, or bounds q·π
static DefiniteResult definite_exp(Expr* f, const Bounds* ab) {
    DefiniteResult r = { NULL, 0, NULL, NULL };
    ExpSumShape sh;
    const char* why = NULL;
    if (!expsum_shape(f, &sh, &why)) {
        if (why == EXPSUM_SUM_DIVISOR) return definite_tower(f, ab);
        r.error = "definite integrals need a rational function, or exp, sin and cos without ln, for now";
        return r;
    }
    Span s = span_of(ab);
    int dir = end_cmp(&ab->a, &ab->b);
    r.error = exp_obstacle(&sh, &s);
    if (!r.error && (dir == 0 || (sh.pure && poly_is_zero(sh.c0_num)))) r.exact = create_const(0);
    else if (!r.error) {
        double scale = 0;
        r.approx = (dir < 0 ? 1 : -1) * numeric_value_expr(f, &s, &scale);
        Expr* F = exp_antiderivative(&r, f, &sh);
        if (F) exact_or_note(&r, &(Problem){ f, ab, &s, scale, NULL, NULL }, F);
    }
    rat_clear(&s.lo); rat_clear(&s.hi);
    expsum_shape_free(&sh);
    return r;
}

DefiniteResult integrate_definite(Expr* f, const Bounds* ab) {
    if (expr_has_exp(f) || expr_has_trig(f) || expr_has_ln(f) || ab->a.pi || ab->b.pi) return definite_exp(f, ab);
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
        r.error = DIVERGES_AT_INF;
    else if (dir == 0 || poly_is_zero(fn.num)) r.exact = create_const(0);
    else if (vanishes_on(fn.den, &s)) r.error = "the integrand has a pole in [a, b], so the integral diverges";
    else {
        r.approx = (dir < 0 ? 1 : -1) * numeric_value(&fn, &s);
        Expr* F = rational_antiderivative(&r, f);
        if (F) exact_or_note(&r, &(Problem){ f, ab, &s, 0, NULL, NULL }, F);
    }
    rat_clear(&s.lo); rat_clear(&s.hi);
    free_polynomial(one);
    free_expr(c);
    return r;
}

// --- bounds from text ---

// inf, +inf, oo, +oo → +1; -inf, -oo → -1; else 0
static int infinity_sign(const char* text) {
    int sign = 1;
    if (*text == '+' || *text == '-') sign = *text++ == '-' ? -1 : 1;
    return strcmp(text, "inf") == 0 || strcmp(text, "oo") == 0 ? sign : 0;
}

// a constant expression ("3/4", "-2") evaluated exactly into *out
static bool parse_constant(Rational* out, const char* text) {
    Expr* e = parse_and_simplify(text);
    bool ok = e && e->type == EXPR_POLY && polynomial_degree(e->poly) <= 0;
    if (ok) poly_coeff(out, e->poly, 0);
    if (e) free_expr(e);
    return ok;
}

// q·π written c·pi/d with optional c ("", "-", "3", "3*") and "/d": pi, -pi/2, 3pi/4, 2*pi
static bool parse_pi_bound(Rational* q, const char* text) {
    const char* p = strstr(text, "pi");
    if (!p) return false;
    size_t n = (size_t)(p - text);
    char* pre = malloc(n + 1);
    memcpy(pre, text, n);
    pre[n] = '\0';
    if (n && pre[n - 1] == '*') pre[--n] = '\0';
    bool ok = true;
    if (n == 0 || strcmp(pre, "+") == 0) rat_set_si(q, 1, 1);
    else if (strcmp(pre, "-") == 0) rat_set_si(q, -1, 1);
    else ok = parse_constant(q, pre);
    free(pre);
    const char* post = p + 2;
    Rational d;
    rat_init(&d);
    if (ok && *post) ok = *post == '/' && parse_constant(&d, post + 1) && !rat_is_zero(&d);
    if (ok && *post) rat_div(q, q, &d);
    rat_clear(&d);
    return ok;
}

bool parse_endpoint(Rational* out, Endpoint* end, const char* text) {
    end->inf = infinity_sign(text);
    end->pi = false;
    if (end->inf) return true;
    if (strstr(text, "pi")) {
        if (!parse_pi_bound(out, text)) return false;
        end->pi = !rat_is_zero(out);                       // 0·π = 0
        return true;
    }
    return parse_constant(out, text);
}

// the integrand: the raw tree with exp, sin, cos or ln (definite.c checks it), else simplified
Expr* parse_integrand(const char* text) {
    Expr* raw = parse_input(text);
    if (raw && (expr_has_exp(raw) || expr_has_ln(raw) || expr_has_trig(raw))) return raw;
    if (raw) free_expr(raw);
    return parse_and_simplify(text);
}

void definite_free(DefiniteResult* r) {
    if (r->exact) free_expr(r->exact);
    r->exact = NULL;
}
