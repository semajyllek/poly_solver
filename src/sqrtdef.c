// Definite integrals with one square root: ∫_a^b (a(x) + b(x)·sqrt(U)) on rational or infinite bounds.
// The range is checked exactly (U >= 0 on [a, b]: no root of U strictly inside, U > 0 at the
// middle; no pole inside; at each end an exact order of a + b·y above -1, cancellation included),
// the antiderivative is evaluated at both ends (algebraic pieces together, and the arguments of
// ln and atan, through Q(x)(y) with limits where a denominator vanishes), atan jumps inside
// corrected, folded exactly, and compared with tanh-sinh quadrature, which never samples an
// endpoint (1/sqrt(1 - x^2) at 1). An infinite end is mapped to u = 0 by x = 1/u.
#include "algsqrt.h"
#include "definite.h"
#include "constfold.h"
#include "roots.h"
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// --- exact checks on [lo, hi] ---

// the square-free part of p
static Polynomial* square_free(const Polynomial* p) {
    Polynomial *dp = poly_derivative(p), *g = poly_gcd(p, dp), *s = poly_divide(p, g);
    free_polynomial(dp); free_polynomial(g);
    return s;
}

typedef struct {
    const Rational *lo, *hi;
    bool closed;            // [lo, hi] or (lo, hi)
} Range;

static bool in_range(const Rational* r, const Range* R) {
    int a = rat_cmp(r, R->lo), b = rat_cmp(r, R->hi);
    return R->closed ? a >= 0 && b <= 0 : a > 0 && b < 0;
}

// a real root of p in the range, decided exactly (isolating intervals refined until each is
// clear of the range or inside it; a root at an end is found by evaluation)
static bool root_in(const Polynomial* p, const Range* R) {
    if (polynomial_degree(p) <= 0) return false;
    Polynomial* s = square_free(p);
    RootIntervalList L = isolate_real_roots(s);
    Rational v, w;
    rat_init(&v); rat_init(&w);
    bool found = false;
    for (const Rational* e = R->lo; e && !found; e = e == R->lo ? R->hi : NULL) {   // roots at the ends
        poly_eval_rat(&v, s, e);
        found = R->closed && rat_is_zero(&v);
    }
    for (int i = 0; i < L.count && !found; i++) {
        RootInterval* I = &L.items[i];
        for (int step = 0; step < 400; step++) {
            if (rat_cmp(&I->lo, &I->hi) == 0) { found = in_range(&I->lo, R); break; }
            if (rat_cmp(&I->hi, R->lo) <= 0 || rat_cmp(&I->lo, R->hi) >= 0) break;       // clear of the range
            poly_eval_rat(&v, s, R->lo); bool at_lo = rat_is_zero(&v) && rat_cmp(&I->lo, R->lo) < 0 && rat_cmp(R->lo, &I->hi) < 0;
            poly_eval_rat(&v, s, R->hi); bool at_hi = rat_is_zero(&v) && rat_cmp(&I->lo, R->hi) < 0 && rat_cmp(R->hi, &I->hi) < 0;
            if (at_lo || at_hi) { found = R->closed; break; }    // the interval's one root is that end
            if (rat_cmp(&I->lo, R->lo) >= 0 && rat_cmp(&I->hi, R->hi) <= 0) { found = true; break; }
            rat_sub(&w, &I->hi, &I->lo);
            mpq_div_2exp(w.val, w.val, 1);
            refine_root(I, s, &w);
        }
    }
    rat_clear(&v); rat_clear(&w);
    root_interval_list_free(&L);
    free_polynomial(s);
    return found;
}

// U >= 0 on [lo, hi]: no root strictly inside and U > 0 at the middle
static bool real_on(const Polynomial* U, const Rational* lo, const Rational* hi) {
    if (!U) return true;
    Rational mid, v;
    rat_init(&mid); rat_init(&v);
    rat_add(&mid, lo, hi);
    mpq_div_2exp(mid.val, mid.val, 1);
    poly_eval_rat(&v, U, &mid);
    bool ok = mpq_sgn(v.val) > 0 && !root_in(U, &(Range){ lo, hi, false });
    rat_clear(&mid); rat_clear(&v);
    return ok;
}

// an²·bd² - bn²·U·ad²: the numerator of a² - b²·U
static Polynomial* conj_numerator(const SqrtShape* s) {
    Polynomial *p = poly_multiply(s->num_a, s->den_b), *q = poly_multiply(s->num_b, s->den_a);
    Polynomial *p2 = poly_multiply(p, p), *q2 = poly_multiply(q, q), *q2u = poly_multiply(q2, s->U), *N = poly_subtract(p2, q2u);
    free_polynomial(p); free_polynomial(q); free_polynomial(p2); free_polynomial(q2); free_polynomial(q2u);
    return N;
}

// the multiplicity of e as a root of p (p ≠ 0), and p/(x - e)^k at e in *v
static int order_at(const Polynomial* p, const Rational* e, Rational* v) {
    Rational one = rat_one_val(), m;
    rat_init(&m);
    rat_neg(&m, e);
    Polynomial *lin = create_polynomial(), *q = poly_copy(p);
    add_term(lin, &one, 1); add_term(lin, &m, 0);
    int k = 0;
    for (poly_eval_rat(v, q, e); rat_is_zero(v); poly_eval_rat(v, q, e), k++) {
        Polynomial* d = poly_divide(q, lin);
        free_polynomial(q);
        q = d;
    }
    free_polynomial(lin); free_polynomial(q);
    rat_clear(&one); rat_clear(&m);
    return k;
}

// twice the order of f = a + b·y at e (f ~ c·(x - e)^(k/2)), exactly: the smaller of the orders of
// a and b·y, or, when their leading terms cancel (a ~ -b·sqrt(U(e))), ord(a² - b²·U) - ord(a);
// INT_MAX for f = 0
static int twice_order(const SqrtShape* s, const Rational* e) {
    Rational va, vb, vu, t, w;
    rat_init(&va); rat_init(&vb); rat_init(&vu); rat_init(&t); rat_init(&w);
    bool has_a = !poly_is_zero(s->num_a), has_b = s->U && !poly_is_zero(s->num_b);
    int oa = has_a ? 2 * (order_at(s->num_a, e, &va) - order_at(s->den_a, e, &t)) : INT_MAX;
    if (has_a) rat_div(&va, &va, &t);
    int ou = has_b ? order_at(s->U, e, &vu) : 0, ob = INT_MAX;
    if (has_b) { ob = 2 * (order_at(s->num_b, e, &vb) - order_at(s->den_b, e, &t)) + ou; rat_div(&vb, &vb, &t); }
    int k = oa < ob ? oa : ob;
    if (has_a && has_b && oa == ob && ou == 0 && rat_is_negative(&va) != rat_is_negative(&vb)) {
        rat_mul(&t, &va, &va);                                // α² = β²·U(e): the leading terms cancel
        rat_mul(&w, &vb, &vb);
        rat_mul(&w, &w, &vu);
        if (rat_eq(&t, &w)) {
            Polynomial* N = conj_numerator(s);
            k = 2 * (order_at(N, e, &t) - 2 * order_at(s->den_a, e, &w) - 2 * order_at(s->den_b, e, &w)) - oa;
            free_polynomial(N);
        }
    }
    rat_clear(&va); rat_clear(&vb); rat_clear(&vu); rat_clear(&t); rat_clear(&w);
    return k;
}

// --- the antiderivative at a point ---

// e depends on x
static bool has_x(const Expr* e) {
    switch (e->type) {
        case EXPR_VAR: return true;
        case EXPR_POLY: return polynomial_degree(e->poly) > 0;
        case EXPR_RATIONAL_FN: return polynomial_degree(e->ratfn.num) > 0 || polynomial_degree(e->ratfn.den) > 0;
        case EXPR_BINOP: return has_x(e->binop.left) || has_x(e->binop.right);
        case EXPR_NEG: return has_x(e->operand);
        case EXPR_LN: return has_x(e->ln_arg);
        case EXPR_ATAN: case EXPR_SQRT: case EXPR_EXP: case EXPR_SIN: case EXPR_COS: return has_x(e->fn_arg);
        case EXPR_ROOTSUM: return true;
        default: return false;
    }
}

typedef struct {
    const Polynomial* U;
    const Rational* x0;
    int inward;             // +1 at the lower end, -1 at the upper: the side the range lies on
} At;

// e with x = x0 everywhere, inside square roots too (expr_at leaves sqrt(U(x)) alone); NULL at a
// zero denominator of a rational function
static Expr* substitute(const Expr* e, const Rational* x0) {
    Rational v, d;
    rat_init(&v); rat_init(&d);
    Expr* r = NULL;
    switch (e->type) {
        case EXPR_VAR: r = create_const_rat(x0); break;
        case EXPR_POLY: poly_eval_rat(&v, e->poly, x0); r = create_const_rat(&v); break;
        case EXPR_RATIONAL_FN:
            poly_eval_rat(&d, e->ratfn.den, x0);
            if (!rat_is_zero(&d)) { poly_eval_rat(&v, e->ratfn.num, x0); rat_div(&v, &v, &d); r = create_const_rat(&v); }
            break;
        case EXPR_BINOP: {
            Expr *l = substitute(e->binop.left, x0), *rr = l ? substitute(e->binop.right, x0) : NULL;
            bool ok = true;
            if (rr && e->binop.op == OP_DIV && expr_eval(rr, 0, &ok) == 0) { free_expr(rr); rr = NULL; }  // 0/0 or c/0: no value here
            if (l && rr) r = e->binop.op == OP_ADD ? create_add(l, rr) : e->binop.op == OP_SUB ? create_sub(l, rr)
                           : e->binop.op == OP_MUL ? create_mul(l, rr) : e->binop.op == OP_DIV ? create_div(l, rr) : create_pow(l, rr);
            else { free_expr(l); free_expr(rr); }
            break;
        }
        case EXPR_NEG: { Expr* o = substitute(e->operand, x0); r = o ? create_neg(o) : NULL; break; }
        case EXPR_SQRT: { Expr* o = substitute(e->fn_arg, x0); r = o ? create_sqrt(o) : NULL; break; }
        case EXPR_LN: { Expr* o = substitute(e->ln_arg, x0); r = o ? create_ln(o) : NULL; break; }
        case EXPR_ATAN: { Expr* o = substitute(e->fn_arg, x0); r = o ? create_atan(o) : NULL; break; }
        default: r = expr_copy(e);                              // constants
    }
    rat_clear(&v); rat_clear(&d);
    return r;
}

static Expr* value_of(const Expr* arg, const At* p);

// a product/quotient chain split into its x-free factors (*k) and the rest (*g), each NULL for none
typedef struct {
    Expr *k, *g;
} Factors;

static void push_factor(Expr** acc, Expr* f, bool inv) {
    if (!*acc) *acc = inv ? create_div(create_const(1), f) : f;
    else *acc = inv ? create_div(*acc, f) : create_mul(*acc, f);
}

static void split_factors(const Expr* e, bool inv, Factors* F) {
    if (e->type == EXPR_BINOP && (e->binop.op == OP_MUL || e->binop.op == OP_DIV)) {
        split_factors(e->binop.left, inv, F);
        split_factors(e->binop.right, e->binop.op == OP_DIV ? !inv : inv, F);
        return;
    }
    if (e->type == EXPR_NEG) { push_factor(&F->k, create_const(-1), false); split_factors(e->operand, inv, F); return; }
    push_factor(has_x(e) ? &F->g : &F->k, expr_copy(e), inv);
}

// e·w with w distributed over e's sums (so each summand's limit is taken with it)
static Expr* times_distributed(const Expr* e, const Polynomial* w) {
    if (e->type == EXPR_NEG) return create_neg(times_distributed(e->operand, w));
    if (e->type == EXPR_BINOP && (e->binop.op == OP_ADD || e->binop.op == OP_SUB)) {
        Expr *l = times_distributed(e->binop.left, w), *r = times_distributed(e->binop.right, w);
        return e->binop.op == OP_ADD ? create_add(l, r) : create_sub(l, r);
    }
    return create_mul(expr_copy(e), create_poly_expr(poly_copy(w)));
}

// A/B at x0 when both are infinite there: (A·(x - x0)^k)/(B·(x - x0)^k) for the least k <= 4 that
// makes both finite with B nonzero; NULL otherwise
static Expr* scaled_quotient(const Expr* arg, const At* p) {
    Rational one = rat_one_val(), m;
    rat_init(&m);
    rat_neg(&m, p->x0);
    Expr* v = NULL;
    for (int k = 1; k <= 4 && !v; k++) {
        Polynomial *lin = create_polynomial(), *w;
        add_term(lin, &one, 1); add_term(lin, &m, 0);
        w = poly_pow(lin, (unsigned)k);
        Expr *A = times_distributed(arg->binop.left, w), *B = times_distributed(arg->binop.right, w);
        free_polynomial(w);
        Expr *a = value_of(A, p), *b = a ? value_of(B, p) : NULL;
        bool ok = true;
        if (b && expr_eval(b, 0, &ok) != 0 && ok) { Expr* q = create_div(a, b); v = expr_at(q, p->x0); free_expr(q); a = b = NULL; }
        free_expr(a); free_expr(b); free_expr(A); free_expr(B);
        free_polynomial(lin);
    }
    rat_clear(&one); rat_clear(&m);
    return v;
}

// arg's value at x0: through Q(x)(y) as a whole (with limits), else its two sides apart (a sum or
// quotient mixing sqrt(3) with y, each side's limit taken in Q(x)(y)) and folded by expr_at, which
// multiplies constant roots into one (sqrt(3)·sqrt(U(x0))); NULL when there is no finite value
static Expr* value_of(const Expr* arg, const At* p) {
    if (!has_x(arg)) return expr_copy(arg);
    if (arg->type == EXPR_NEG) { Expr* v = value_of(arg->operand, p); return v ? create_neg(v) : NULL; }
    Expr* v = sqrt_value_at(arg, p->U, p->x0);
    if (v) return v;
    if (arg->type == EXPR_BINOP && (arg->binop.op == OP_MUL || arg->binop.op == OP_DIV)) {   // c·g: g's limit through Q(x)(y)
        Factors F = { NULL, NULL };
        split_factors(arg, false, &F);
        Expr* g = F.k && F.g ? sqrt_value_at(F.g, p->U, p->x0) : NULL;
        if (g) { Expr* e = create_mul(F.k, g); v = expr_at(e, p->x0); free_expr(e); free_expr(F.g); return v; }
        free_expr(F.k); free_expr(F.g);
    }
    if (arg->type == EXPR_BINOP && arg->binop.op != OP_POW) {
        bool div = arg->binop.op == OP_DIV;
        Expr *l = value_of(arg->binop.left, p), *r = l || div ? value_of(arg->binop.right, p) : NULL;
        if (div && (!l || !r)) { free_expr(l); free_expr(r); return scaled_quotient(arg, p); }   // ∞/∞ or ∞/c
        bool ok = true;
        if (r && div && expr_eval(r, 0, &ok) == 0) { free_expr(r); r = NULL; }   // c/0: no value
        if (!r) { free_expr(l); return NULL; }
        BinOp op = arg->binop.op;
        Expr* e = op == OP_ADD ? create_add(l, r) : op == OP_SUB ? create_sub(l, r) : op == OP_MUL ? create_mul(l, r) : create_div(l, r);
        v = expr_at(e, p->x0);
        free_expr(e);
        return v;
    }
    Expr* s = substitute(arg, p->x0);
    if (!s) return NULL;
    v = expr_at(s, p->x0);
    free_expr(s);
    return v;
}

// atan(arg) at an end where |arg| → ∞: ±π/2 with the sign of arg just inside the range
static Expr* atan_at_infinity(const Expr* arg, const At* p) {
    double x = rat_to_double(p->x0), s = 1 + fabs(x);
    bool ok = true;
    double v1 = expr_eval(arg, x + p->inward * 1e-6 * s, &ok), v = expr_eval(arg, x + p->inward * 1e-10 * s, &ok);
    if (!ok || !isfinite(v) || v * v1 <= 0 || fabs(v) < 10 * fabs(v1) || fabs(v) < 1e3) return NULL; // no clear approach to ±∞
    Rational half;
    rat_init(&half);
    mpq_set_si(half.val, v > 0 ? 1 : -1, 2);
    Expr* e = create_mul(create_const_rat(&half), create_pi());
    rat_clear(&half);
    return e;
}

// one term of F at x0, constant factors kept, ln and atan of their arguments' values, and the
// algebraic rest through Q(x)(y); NULL when a value is infinite (a log at a zero) or unknown
static Expr* term_at(const Expr* t, const At* p) {
    if (!has_x(t)) return expr_copy(t);
    if (t->type == EXPR_NEG) { Expr* v = term_at(t->operand, p); return v ? create_neg(v) : NULL; }
    if (t->type == EXPR_BINOP && t->binop.op == OP_MUL && (!has_x(t->binop.left) || !has_x(t->binop.right))) {
        bool left = !has_x(t->binop.left);
        Expr* v = term_at(left ? t->binop.right : t->binop.left, p);
        if (!v) return NULL;
        return left ? create_mul(expr_copy(t->binop.left), v) : create_mul(v, expr_copy(t->binop.right));
    }
    if (t->type == EXPR_LN || t->type == EXPR_ATAN) {
        const Expr* arg = t->type == EXPR_LN ? t->ln_arg : t->fn_arg;
        Expr* v = value_of(arg, p);
        if (!v) return t->type == EXPR_ATAN ? atan_at_infinity(arg, p) : NULL;
        return t->type == EXPR_LN ? create_ln(v) : create_atan(v);
    }
    return value_of(t, p);
}

// --- atan jumps inside the range ---

// the atan factor of a product term (through signs and constant factors), if any
static const Expr* atan_factor(const Expr* t) {
    if (t->type == EXPR_ATAN) return t;
    if (t->type == EXPR_NEG) return atan_factor(t->operand);
    if (t->type != EXPR_BINOP || t->binop.op != OP_MUL) return NULL;
    const Expr* l = atan_factor(t->binop.left);
    return l ? l : atan_factor(t->binop.right);
}

// the nonconstant polynomial factor of a product (through signs), if any
static const Polynomial* poly_factor(const Expr* e) {
    if (e->type == EXPR_POLY) return polynomial_degree(e->poly) > 0 ? e->poly : NULL;
    if (e->type == EXPR_NEG) return poly_factor(e->operand);
    if (e->type != EXPR_BINOP || e->binop.op != OP_MUL) return NULL;
    const Polynomial* l = poly_factor(e->binop.left);
    return l ? l : poly_factor(e->binop.right);
}

// the polynomial W in an atan argument c·N/(W·y), if any: its zeros are the argument's poles
static const Polynomial* atan_pole_poly(const Expr* arg) {
    if (arg->type == EXPR_NEG) return atan_pole_poly(arg->operand);
    if (arg->type != EXPR_BINOP) return NULL;
    if (arg->binop.op == OP_DIV) return poly_factor(arg->binop.right);
    if (arg->binop.op != OP_MUL) return NULL;
    const Polynomial* l = has_x(arg->binop.left) ? atan_pole_poly(arg->binop.left) : NULL;
    return l ? l : has_x(arg->binop.right) ? atan_pole_poly(arg->binop.right) : NULL;
}

static int sign_at(const Expr* arg, double x) {
    bool ok = true;
    double v = expr_eval(arg, x, &ok);
    return ok ? (v > 0) - (v < 0) : 0;
}

// Σ over the zeros of W strictly inside (lo, hi) of the jump of atan(arg) there, in units of π
// (±∞ on either side: ±π/2 → ∓π/2); each zero isolated exactly, the signs read just beside it
static int atan_jumps(const Expr* arg, const Rational* lo, const Rational* hi_or_inf) {
    const Polynomial* W = atan_pole_poly(arg);
    if (!W) return 0;
    Polynomial* sw = square_free(W);
    Rational B;                                                  // ∞: Cauchy's bound 1 + max |w_i/w_n| on the roots
    rat_init(&B);
    if (!hi_or_inf) {
        Rational l, c;
        rat_init(&l); rat_init(&c);
        poly_coeff(&l, sw, polynomial_degree(sw));
        for (const Term* t = sw->head; t; t = t->next) {
            rat_div(&c, &t->coeff, &l);
            rat_abs(&c, &c);
            if (rat_cmp(&c, &B) > 0) rat_set(&B, &c);
        }
        mpq_set_ui(c.val, 1, 1);
        rat_add(&B, &B, &c);
        rat_clear(&l); rat_clear(&c);
    }
    const Rational* hi = hi_or_inf ? hi_or_inf : &B;
    RootIntervalList L = isolate_real_roots(sw);
    Rational w;
    rat_init(&w);
    int jumps = 0;
    for (int i = 0; i < L.count; i++) {
        RootInterval* I = &L.items[i];
        for (int step = 0; step < 200; step++) {
            bool exact = rat_cmp(&I->lo, &I->hi) == 0;
            if (rat_cmp(&I->hi, lo) <= 0 || rat_cmp(&I->lo, hi) >= 0) break;            // outside (or at an end)
            if (exact || (rat_cmp(&I->lo, lo) > 0 && rat_cmp(&I->hi, hi) < 0 && step >= 40)) {
                if (exact && (rat_cmp(&I->lo, lo) <= 0 || rat_cmp(&I->lo, hi) >= 0)) break;
                double a = rat_to_double(&I->lo), b = rat_to_double(&I->hi);
                if (exact) { a -= 1e-9 * (1 + fabs(a)); b += 1e-9 * (1 + fabs(b)); }
                jumps += (sign_at(arg, b) - sign_at(arg, a)) / 2;
                break;
            }
            rat_sub(&w, &I->hi, &I->lo);
            mpq_div_2exp(w.val, w.val, 1);
            refine_root(I, sw, &w);
        }
    }
    rat_clear(&w); rat_clear(&B);
    root_interval_list_free(&L);
    free_polynomial(sw);
    return jumps;
}

// the ln factor of a product term (through signs and constant factors), if any
static const Expr* ln_factor(const Expr* t) {
    if (t->type == EXPR_LN) return t;
    if (t->type == EXPR_NEG) return ln_factor(t->operand);
    if (t->type != EXPR_BINOP || t->binop.op != OP_MUL) return NULL;
    const Expr* l = ln_factor(t->binop.left);
    return l ? l : ln_factor(t->binop.right);
}

// copy of the product t with its ln factor f replaced by g (taken)
static Expr* replaced_ln(const Expr* t, const Expr* f, Expr* g) {
    if (t == f) return g;
    if (t->type == EXPR_NEG) return create_neg(replaced_ln(t->operand, f, g));
    bool left = ln_factor(t->binop.left) == f;
    return left ? create_mul(replaced_ln(t->binop.left, f, g), expr_copy(t->binop.right))
                : create_mul(expr_copy(t->binop.left), replaced_ln(t->binop.right, f, g));
}

// copy of the product t with its factor f replaced by k·π
static Expr* with_pi(const Expr* t, const Expr* f, long k) {
    if (t == f) return create_mul(create_const(k), create_pi());
    if (t->type == EXPR_NEG) return create_neg(with_pi(t->operand, f, k));
    bool left = atan_factor(t->binop.left) == f;
    return left ? create_mul(with_pi(t->binop.left, f, k), expr_copy(t->binop.right))
                : create_mul(expr_copy(t->binop.left), with_pi(t->binop.right, f, k));
}

// the summands of F
static void summands(const Expr* F, int sign, Expr*** out, int* n) {
    if (F->type == EXPR_BINOP && (F->binop.op == OP_ADD || F->binop.op == OP_SUB)) {
        summands(F->binop.left, sign, out, n);
        summands(F->binop.right, F->binop.op == OP_ADD ? sign : -sign, out, n);
        return;
    }
    *out = realloc(*out, (size_t)(*n + 1) * sizeof(Expr*));
    (*out)[(*n)++] = sign > 0 ? expr_copy(F) : create_neg(expr_copy(F));
}

// an ln or atan inside e
static bool has_transcendental(const Expr* e) {
    switch (e->type) {
        case EXPR_LN: case EXPR_ATAN: case EXPR_ROOTSUM: return true;
        case EXPR_BINOP: return has_transcendental(e->binop.left) || has_transcendental(e->binop.right);
        case EXPR_NEG: return has_transcendental(e->operand);
        default: return false;
    }
}

// the summands of F, the algebraic ones merged into one (its pieces may each be infinite at an
// end, as (x² + 2/3)/x³ and (2/3·x² + 2/3)·y/x³ at 0, while their sum has a limit)
static void grouped_summands(const Expr* F, Expr*** T, int* n) {
    Expr** all = NULL;
    int m = 0;
    summands(F, 1, &all, &m);
    Expr* alg = NULL;
    *T = malloc((size_t)(m + 1) * sizeof(Expr*));
    *n = 0;
    for (int i = 0; i < m; i++) {
        if (has_transcendental(all[i]) || !has_x(all[i])) (*T)[(*n)++] = all[i];
        else alg = alg ? create_add(alg, all[i]) : all[i];
    }
    if (alg) (*T)[(*n)++] = alg;
    free(all);
}

// owned constant terms, summed at the end
typedef struct {
    Expr** t;
    int n;
} Parts;

static void parts_push(Parts* P, Expr* e) {
    P->t = realloc(P->t, (size_t)(P->n + 1) * sizeof(Expr*));
    P->t[P->n++] = e;
}

static void parts_free(Parts* P) {
    for (int i = 0; i < P->n; i++) free_expr(P->t[i]);
    free(P->t);
}

// sign·F(x0) for F at p, into P
typedef struct {
    const Expr* F;
    At p;
    int sign;
} Value;

// sign·F(x0) into P, term by term (the algebraic summands as one); false when a term has no value
static bool add_values(const Value* V, Parts* P) {
    const Expr* F = V->F;
    const At* p = &V->p;
    int sign = V->sign;
    Expr** T = NULL;
    int n = 0;
    grouped_summands(F, &T, &n);
    bool ok = !expr_has_rootsum(F);
    for (int i = 0; i < n && ok; i++) {
        Expr* v = term_at(T[i], p);
        Expr* f = v ? expr_at(v, p->x0) : NULL;
        free_expr(v);
        ok = f != NULL;
        if (ok) parts_push(P, sign > 0 ? f : create_neg(f));
    }
    for (int i = 0; i < n; i++) free_expr(T[i]);
    free(T);
    return ok;
}

// (lo, hi), hi NULL for ∞, and the sign the jumps enter with
typedef struct {
    const Rational *lo, *hi;
    int sign;
} Jumps;

// sign·(the atan jumps of F's terms on (lo, hi)) in multiples of π into P, so that F minus its
// jumps is continuous there
static void add_jumps(const Expr* F, const Jumps* J, Parts* P) {
    Expr** T = NULL;
    int n = 0;
    summands(F, 1, &T, &n);
    for (int i = 0; i < n; i++) {
        const Expr* at = atan_factor(T[i]);
        int k = at ? atan_jumps(at->fn_arg, J->lo, J->hi) : 0;
        if (k) parts_push(P, with_pi(T[i], at, J->sign * k));
        free_expr(T[i]);
    }
    free(T);
}

// the folded sum of P (taken), or NULL when ok is false
static Expr* parts_sum(Parts* P, bool ok) {
    if (!ok) { parts_free(P); return NULL; }
    Expr* r = const_fold_sum(P->t, P->n);
    free(P->t);
    return r;
}

// F(b) - F(a), each term evaluated and folded, minus the atan jumps crossed; NULL when a term has
// no value at an end
static Expr* exact_difference(const Expr* F, const Polynomial* U, const Rational* const* ends) {
    bool up = rat_cmp(ends[0], ends[1]) < 0;
    const Rational *lo = up ? ends[0] : ends[1], *hi = up ? ends[1] : ends[0];
    Parts P = { NULL, 0 };
    bool ok = add_values(&(Value){ F, { U, ends[1], up ? -1 : 1 }, 1 }, &P) && add_values(&(Value){ F, { U, ends[0], up ? 1 : -1 }, -1 }, &P);
    if (ok) add_jumps(F, &(Jumps){ lo, hi, up ? -1 : 1 }, &P);
    return parts_sum(&P, ok);
}

// --- tanh-sinh quadrature, from each end ---

// f at e + s·u for u >= 0 into the range (s = +1 at the lower end, -1 at the upper): the
// polynomials of a + b·sqrt(U) shifted exactly to e, so u near 0 (a singular end) loses nothing
typedef struct {
    DPoly an, ad, bn, bd, U;
    DPoly N, D, P, Q;   // the conjugate form u^eN·N/(D·(u^eP·P - u^eQ·Q·y)) with powers of u divided out:
    int eN, eP, eQ;     // N ~ an²·bd² - bn²·U·ad², D ~ ad·bd, P ~ an·bd, Q ~ bn·ad
    bool root;
} Near;

// p(e + s·u) as a polynomial in u, exactly
static Polynomial* shift_exact(const Polynomial* p, const Rational* e, int s) {
    Rational sr = rat_from_int(s);
    Polynomial *lin = create_polynomial(), *r = create_polynomial();
    add_term(lin, &sr, 1); add_term(lin, e, 0);
    Rational c;
    rat_init(&c);
    for (int k = polynomial_degree(p); k >= 0; k--) {         // Horner
        poly_coeff(&c, p, k);
        Polynomial* t = poly_multiply(r, lin);
        add_term(t, &c, 0);
        free_polynomial(r);
        r = t;
    }
    free_polynomial(lin);
    rat_clear(&c); rat_clear(&sr);
    return r;
}

// the same in doubles
static DPoly shifted(const Polynomial* p, const Rational* e, int s) {
    Polynomial* r = shift_exact(p, e, s);
    DPoly d = dpoly_from(r);
    free_polynomial(r);
    return d;
}

// an end e and the side the range lies on (+1 above, -1 below)
typedef struct {
    const Rational* e;
    int side;
} End;

// p(e + side·u)/u^k in doubles, k the lowest power of u in it (0 for p = 0)
static DPoly shifted_stripped(const Polynomial* p, const End* at, int* k) {
    Polynomial* r = shift_exact(p, at->e, at->side);
    *k = INT_MAX;
    for (const Term* t = r->head; t; t = t->next) if (t->exponent < *k) *k = t->exponent;
    if (*k == INT_MAX) *k = 0;
    Polynomial* q = create_polynomial();
    for (const Term* t = r->head; t; t = t->next) add_term(q, &t->coeff, t->exponent - *k);
    DPoly d = dpoly_from(q);
    free_polynomial(r); free_polynomial(q);
    return d;
}

static Near near_end(const SqrtShape* s, const Rational* e, int side) {
    Near n = { shifted(s->num_a, e, side), shifted(s->den_a, e, side), shifted(s->num_b, e, side), shifted(s->den_b, e, side),
               { NULL, -1 }, { NULL, -1 }, { NULL, -1 }, { NULL, -1 }, { NULL, -1 }, 0, 0, 0, s->U != NULL };
    if (!s->U) return n;
    Polynomial *N = conj_numerator(s), *D = poly_multiply(s->den_a, s->den_b);
    Polynomial *P = poly_multiply(s->num_a, s->den_b), *Q = poly_multiply(s->num_b, s->den_a);
    int kN, kD, kP, kQ;
    n.U = shifted(s->U, e, side);
    End at = { e, side };
    n.N = shifted_stripped(N, &at, &kN); n.D = shifted_stripped(D, &at, &kD);
    n.P = shifted_stripped(P, &at, &kP); n.Q = shifted_stripped(Q, &at, &kQ);
    int m = poly_is_zero(P) ? kQ : poly_is_zero(Q) ? kP : (kP < kQ ? kP : kQ);
    n.eN = kN - kD - m; n.eP = kP - m; n.eQ = kQ - m;
    free_polynomial(N); free_polynomial(D); free_polynomial(P); free_polynomial(Q);
    return n;
}

static void near_free(Near* n) {
    dpoly_free(&n->an); dpoly_free(&n->ad); dpoly_free(&n->bn); dpoly_free(&n->bd);
    if (n->root) { dpoly_free(&n->U); dpoly_free(&n->N); dpoly_free(&n->D); dpoly_free(&n->P); dpoly_free(&n->Q); }
}

// a + b·y at u, through the conjugate where the two nearly cancel or overflow ((sqrt(x² + 1) - x)² at x = 1/u)
static double near_at(const Near* n, double u) {
    double v = dpoly_eval(&n->an, u) / dpoly_eval(&n->ad, u);
    if (!n->root) return v;
    double y = sqrt(fmax(dpoly_eval(&n->U, u), 0));            // U >= 0 here
    double w = dpoly_eval(&n->bn, u) / dpoly_eval(&n->bd, u) * y;
    if (isfinite(v + w) && fabs(v + w) >= 1e-3 * (fabs(v) + fabs(w))) return v + w;
    double den = dpoly_eval(&n->D, u) * (pow(u, n->eP) * dpoly_eval(&n->P, u) - pow(u, n->eQ) * dpoly_eval(&n->Q, u) * y);
    return pow(u, n->eN) * dpoly_eval(&n->N, u) / den;
}

// ∫_0^H g(u) du by tanh-sinh, u = H/(e^(-2v) + 1) with v = π/2·sinh t: the distance to u = 0
// (the end) computed directly, halving the step until two estimates agree
static double tanh_sinh(const Near* n, double H) {
    double prev = NAN;
    for (int level = 0; level < 12; level++) {
        double h = ldexp(1, -level), sum = 0;
        for (double t = -4; t <= 4 + 1e-12; t += h) {
            double v = M_PI / 2 * sinh(t), ch = cosh(v), w = H / 2 * M_PI / 2 * cosh(t) / (ch * ch);
            double u = t < 0 ? H / (exp(-2 * v) + 1) : H - H / (exp(2 * v) + 1);
            if (w < 1e-300 || u <= 0) continue;                // weights underflow far out; an end itself
            double g = near_at(n, u);
            if (!isfinite(g)) return NAN;
            sum += w * g;
        }
        double est = sum * h;
        if (!isnan(prev) && fabs(est - prev) <= 1e-13 * (1 + fabs(est))) return est;
        prev = est;
    }
    return prev;
}

// ∫_lo^hi f: the two halves, each from its own end
static double quadrature(const SqrtShape* s, const Rational* lo, const Rational* hi) {
    double H = (rat_to_double(hi) - rat_to_double(lo)) / 2;
    Near a = near_end(s, lo, 1), b = near_end(s, hi, -1);
    double v = tanh_sinh(&a, H) + tanh_sinh(&b, H);
    near_free(&a); near_free(&b);
    return v;
}

// --- entry point ---

#define AGREE_SQRT 1e-8 // relative agreement required between the exact value and the quadrature

static DefiniteResult fail(const char* error) { return (DefiniteResult){ NULL, 0, error, NULL }; }

static const char* const POLE_ERROR = "the integrand has a pole in [a, b], so the integral diverges";

// why ∫_lo^hi f has no value (sqrt(U) not real, or a pole that is not integrable), or NULL
static const char* range_error(const SqrtShape* s, const Rational* lo, const Rational* hi) {
    if (!real_on(s->U, lo, hi)) return "the square root is not real on all of [a, b]";
    if (root_in(s->den_a, &(Range){ lo, hi, false }) || root_in(s->den_b, &(Range){ lo, hi, false })
        || twice_order(s, lo) <= -2 || twice_order(s, hi) <= -2)
        return POLE_ERROR;
    return NULL;
}

// ∫ over finite rational bounds
static DefiniteResult finite_range(const Expr* f, const Rational* a, const Rational* b) {
    const Rational* lo = rat_cmp(a, b) <= 0 ? a : b, *hi = lo == a ? b : a;
    if (rat_cmp(a, b) == 0) return (DefiniteResult){ create_const(0), 0, NULL, NULL };
    SqrtShape s;
    const char* why = NULL;
    if (!sqrt_shape(f, &s, &why)) return fail(why);
    DefiniteResult r = { NULL, 0, range_error(&s, lo, hi), NULL };
    if (!r.error) {
        r.approx = quadrature(&s, lo, hi) * (lo == a ? 1 : -1);
        RischResult F = integrate_sqrt(f);
        Expr* v = F.status == RISCH_OK ? exact_difference(F.result, s.U, (const Rational*[]){ a, b }) : NULL;
        bool ok = true;
        double d = v ? expr_eval(v, 0, &ok) : 0;
        if (F.status != RISCH_OK) r.note = F.reason;
        else if (!v) r.note = "the antiderivative could not be evaluated exactly at a bound";
        else if (ok && fabs(d - r.approx) <= AGREE_SQRT * (1 + fabs(r.approx))) r.exact = v;
        else r.note = "internal check failed: exact value disagrees with quadrature";
        if (v && !r.exact) free_expr(v);
        free_expr(F.result);
    }
    sqrt_shape_free(&s);
    return r;
}

// --- infinite bounds: x = 1/u maps [c, ∞) to (0, 1/c], sqrt(U(1/u)) = sqrt(V(u))/u^k ---

// p(-x)
static Polynomial* mirrored(const Polynomial* p) {
    Polynomial* r = create_polynomial();
    Rational c;
    rat_init(&c);
    for (const Term* t = p->head; t; t = t->next) {
        rat_set(&c, &t->coeff);
        if (t->exponent % 2) rat_neg(&c, &c);
        add_term(r, &c, t->exponent);
    }
    rat_clear(&c);
    return r;
}

// u^n·p(1/u) (n >= deg p)
static Polynomial* reversed(const Polynomial* p, int n) {
    Polynomial* r = create_polynomial();
    for (const Term* t = p->head; t; t = t->next) add_term(r, &t->coeff, n - t->exponent);
    return r;
}

static Polynomial* x_power(int k) {
    Polynomial* p = create_polynomial();
    Rational one = rat_one_val();
    add_term(p, &one, k);
    rat_clear(&one);
    return p;
}

// num(1/u)/(den(1/u)·u^k) as num/den polynomials in u
static Expr* reciprocal_fn(const Polynomial* num, const Polynomial* den, int k) {
    int dn = polynomial_degree(num), dd = polynomial_degree(den);
    Polynomial *n = reversed(num, dn), *d = reversed(den, dd), *un = x_power(dd), *ud = x_power(dn + k);
    Polynomial *N = poly_multiply(n, un), *D = poly_multiply(d, ud);
    free_polynomial(n); free_polynomial(d); free_polynomial(un); free_polynomial(ud);
    return create_rational_fn(N, D);
}

// the shape of f(-x) (in place)
static void mirror_shape(SqrtShape* s) {
    Polynomial** ps[] = { &s->U, &s->num_a, &s->den_a, &s->num_b, &s->den_b };
    for (int i = 0; i < 5; i++)
        if (*ps[i]) { Polynomial* m = mirrored(*ps[i]); free_polynomial(*ps[i]); *ps[i] = m; }
}

// a + b·sqrt(U) as an expression
static Expr* shape_expr(const SqrtShape* s) {
    Expr* a = create_rational_fn(poly_copy(s->num_a), poly_copy(s->den_a));
    if (!s->U) return a;
    return create_add(a, create_mul(create_rational_fn(poly_copy(s->num_b), poly_copy(s->den_b)),
                                    create_sqrt(create_poly_expr(poly_copy(s->U)))));
}

// G(u) = f(1/u)/u² = a(1/u)/u² + b(1/u)/u^(2+k)·sqrt(V(u)), V = u^(2k)·U(1/u), k = ⌈deg U/2⌉
static Expr* tail_integrand(const SqrtShape* s) {
    Expr* A = reciprocal_fn(s->num_a, s->den_a, 2);
    if (!s->U) return A;
    int n = polynomial_degree(s->U), k = (n + 1) / 2;
    Expr* B = reciprocal_fn(s->num_b, s->den_b, 2 + k);
    return create_add(A, create_mul(B, create_sqrt(create_poly_expr(reversed(s->U, 2 * k)))));
}

// e(1/u) with sqrt(p(1/u)) = sqrt(u^(2k)·p(1/u))/u^k (k = ⌈deg p/2⌉); NULL for a RootSum, exp, sin, ...
static Expr* at_reciprocal(const Expr* e) {
    Rational one = rat_one_val();
    Polynomial* unit = poly_constant(&one);
    rat_clear(&one);
    Expr* r = NULL;
    switch (e->type) {
        case EXPR_CONST: case EXPR_TRIG: r = expr_copy(e); break;
        case EXPR_VAR: { Polynomial* u = x_power(1); r = create_rational_fn(poly_copy(unit), u); break; }
        case EXPR_POLY: r = reciprocal_fn(e->poly, unit, 0); break;
        case EXPR_RATIONAL_FN: r = reciprocal_fn(e->ratfn.num, e->ratfn.den, 0); break;
        case EXPR_NEG: { Expr* o = at_reciprocal(e->operand); r = o ? create_neg(o) : NULL; break; }
        case EXPR_LN: { Expr* o = at_reciprocal(e->ln_arg); r = o ? create_ln(o) : NULL; break; }
        case EXPR_ATAN: { Expr* o = at_reciprocal(e->fn_arg); r = o ? create_atan(o) : NULL; break; }
        case EXPR_SQRT:
            if (!has_x(e->fn_arg)) r = expr_copy(e);
            else if (e->fn_arg->type == EXPR_POLY) {
                int k = (polynomial_degree(e->fn_arg->poly) + 1) / 2;
                r = create_div(create_sqrt(create_poly_expr(reversed(e->fn_arg->poly, 2 * k))), create_poly_expr(x_power(k)));
            }
            break;
        case EXPR_BINOP: {
            Expr *l = at_reciprocal(e->binop.left), *rr = l ? at_reciprocal(e->binop.right) : NULL;
            if (!rr) { free_expr(l); break; }
            BinOp op = e->binop.op;
            r = op == OP_ADD ? create_add(l, rr) : op == OP_SUB ? create_sub(l, rr) : op == OP_MUL ? create_mul(l, rr)
              : op == OP_DIV ? create_div(l, rr) : create_pow(l, rr);
            break;
        }
        default: break;
    }
    free_polynomial(unit);
    return r;
}

// the coefficient of a term c·ln|A|, numerically (the term at a point where ln|A| = 1 is not
// available, so the term is divided by its ln's value at x = 1/3); NAN when not of that form
static double ln_coefficient(const Expr* t) {
    const Expr* L = ln_factor(t);
    if (!L) return NAN;
    bool ok = true;
    double v = expr_eval(t, 1.0 / 3, &ok), l = expr_eval(L, 1.0 / 3, &ok);
    return ok && fabs(l) > 1e-9 ? v / l : NAN;
}

// F with each pair c·ln|A| - c·ln|B| merged into c·ln|A/B|: the two may be infinite at ∞ while the
// difference has a limit (-ln|5x² + 2 - 2√2·x·y| + ln|x² + 2|)
static Expr* merged_logs(const Expr* F) {
    Expr** T = NULL;
    int n = 0;
    summands(F, 1, &T, &n);
    Expr* out = NULL;
    for (int i = 0; i < n; i++) {
        if (!T[i]) continue;
        double ci = ln_coefficient(T[i]);
        for (int j = i + 1; j < n && isfinite(ci); j++) {
            double cj = T[j] ? ln_coefficient(T[j]) : NAN;
            if (!isfinite(cj) || fabs(ci + cj) > 1e-12 * fabs(ci)) continue;
            const Expr *A = ln_factor(T[i]), *B = ln_factor(T[j]);   // ci·ln|A| + cj·ln|B| = ci·ln|A/B|
            Expr* q = create_ln(create_div(expr_copy(A->ln_arg), expr_copy(B->ln_arg)));
            Expr* m = replaced_ln(T[i], A, q);
            free_expr(T[i]); free_expr(T[j]);
            T[i] = m; T[j] = NULL;
            break;
        }
        out = expr_append_sum(out, T[i]);
    }
    free(T);
    return out ? out : create_const(0);
}

// ∫_c^∞ f = lim F - F(c) for one antiderivative F, the limit as F(1/u) at u → 0+, minus the atan
// jumps on (c, ∞); NULL with *why when F is not found or a value is not
static Expr* exact_tail(const SqrtShape* s, const Rational* c, const char** why) {
    Expr* f = shape_expr(s);
    RischResult F = integrate_sqrt(f);
    free_expr(f);
    *why = F.status == RISCH_OK ? "the antiderivative could not be evaluated exactly at infinity" : F.reason;
    if (F.status != RISCH_OK) return NULL;
    Expr* M = merged_logs(F.result);
    Expr* Fu = at_reciprocal(M);
    free_expr(M);
    int k = (polynomial_degree(s->U) + 1) / 2;
    Polynomial* V = reversed(s->U, 2 * k);
    Rational zero;
    rat_init(&zero);
    Parts P = { NULL, 0 };
    bool ok = Fu && add_values(&(Value){ Fu, { V, &zero, 1 }, 1 }, &P) && add_values(&(Value){ F.result, { s->U, c, 1 }, -1 }, &P);
    if (ok) add_jumps(F.result, &(Jumps){ c, NULL, -1 }, &P);
    free_expr(Fu); free_expr(F.result);
    free_polynomial(V);
    rat_clear(&zero);
    return parts_sum(&P, ok);
}

// p + q: exact when both are, the error or note of either otherwise (takes both)
static DefiniteResult combine(DefiniteResult p, DefiniteResult q) {
    if (p.error || q.error) {
        DefiniteResult e = p.error ? p : q, o = p.error ? q : p;
        definite_free(&o);
        free_expr(e.exact);
        return (DefiniteResult){ NULL, 0, e.error, NULL };
    }
    DefiniteResult r = { NULL, p.approx + q.approx, NULL, p.note ? p.note : q.note };
    if (p.exact && q.exact) r.exact = const_fold_sum((Expr*[]){ p.exact, q.exact }, 2);
    else { free_expr(p.exact); free_expr(q.exact); }
    return r;
}

// ∫_c^∞ of the shape s, the number from ∫_0^(1/c) G (c > 0) or ∫_c^1 f + ∫_0^1 G, the exact value
// from one antiderivative's limit at ∞ when it agrees; an error when G is not integrable at u = 0
// (f decays no faster than 1/x) or a piece has no value
static DefiniteResult tail(const SqrtShape* s, const Rational* c) {
    Expr* G = tail_integrand(s);
    SqrtShape g;
    const char* why = NULL;
    Rational zero, one, top;
    rat_init(&zero); rat_init(&one); rat_init(&top);
    mpq_set_ui(one.val, 1, 1);
    bool split = mpq_sgn(c->val) <= 0;
    if (split) rat_set(&top, &one);
    else mpq_inv(top.val, c->val);
    DefiniteResult r = { NULL, 0, NULL, NULL };
    if (!sqrt_shape(G, &g, &why)) r.error = why;
    else if (!(r.error = range_error(&g, &zero, &top)) && !(r.error = split ? range_error(s, c, &one) : NULL)) {
        r.approx = quadrature(&g, &zero, &top) + (split ? quadrature(s, c, &one) : 0);
        const char* none = "the antiderivative could not be evaluated exactly at infinity";
        Expr* e = s->U ? exact_tail(s, c, &none) : NULL;
        bool ok = true;
        double d = e ? expr_eval(e, 0, &ok) : NAN;
        if (e && ok && fabs(d - r.approx) <= AGREE_SQRT * (1 + fabs(r.approx))) r.exact = e;
        else {
            free_expr(e);
            r.note = e ? "internal check failed: exact value disagrees with quadrature" : none;
        }
    }
    if (r.error == POLE_ERROR && !why && twice_order(&g, &zero) <= -2) r.error = "the integrand does not decay faster than 1/x, so the integral diverges at infinity";
    if (!why) sqrt_shape_free(&g);
    free_expr(G);
    rat_clear(&zero); rat_clear(&one); rat_clear(&top);
    return r;
}

// -∞ < c < ∞ as -1, 0, 1 for comparison
static int end_rank(const Endpoint* e) { return e->inf; }

// a bound at ±∞: ∫_c^∞ f, ∫_-∞^c f = ∫_-c^∞ f(-x), or both halves at 0; the sign for a > b
static DefiniteResult infinite_range(const Expr* f, const Bounds* ab) {
    bool down = end_rank(&ab->a) > end_rank(&ab->b) || (!ab->a.inf && !ab->b.inf && rat_cmp(ab->a.x, ab->b.x) > 0);
    const Endpoint *lo = down ? &ab->b : &ab->a, *hi = down ? &ab->a : &ab->b;
    if (lo->inf == hi->inf) return (DefiniteResult){ create_const(0), 0, NULL, NULL };   // ∫_∞^∞
    SqrtShape s;
    const char* why = NULL;
    if (!sqrt_shape(f, &s, &why)) return fail(why);
    Rational zero, mc;
    rat_init(&zero); rat_init(&mc);
    DefiniteResult r;
    if (lo->inf && hi->inf) {
        DefiniteResult right = tail(&s, &zero);
        mirror_shape(&s);
        r = combine(right, tail(&s, &zero));
    } else if (hi->inf) r = tail(&s, lo->x);
    else {
        rat_neg(&mc, hi->x);
        mirror_shape(&s);
        r = tail(&s, &mc);
    }
    if (down && !r.error) {
        r.approx = -r.approx;
        if (r.exact) r.exact = const_fold_sum((Expr*[]){ create_neg(r.exact) }, 1);
    }
    sqrt_shape_free(&s);
    rat_clear(&zero); rat_clear(&mc);
    return r;
}

DefiniteResult integrate_definite_sqrt(const Expr* f, const Bounds* ab) {
    if (ab->a.pi || ab->b.pi) return fail("multiples of pi as bounds are not supported with square roots");
    if (ab->a.inf || ab->b.inf) return infinite_range(f, ab);
    return finite_range(f, ab->a.x, ab->b.x);
}
