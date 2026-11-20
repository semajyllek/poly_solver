// Definite integrals with one square root: ∫_a^b (a(x) + b(x)·sqrt(U)) on finite rational bounds.
// The range is checked exactly (U >= 0 on [a, b]: no root of U strictly inside, U > 0 at the
// middle; no pole; an endpoint at a root of U only where the integrand is O((x - e)^(-1/2))),
// the antiderivative is evaluated term by term at both ends (algebraic pieces, and the arguments
// of ln and atan, through Q(x)(y) with limits where a denominator vanishes), folded exactly, and
// compared with tanh-sinh quadrature, which never samples an endpoint (1/sqrt(1 - x^2) at 1).
#include "algsqrt.h"
#include "definite.h"
#include "constfold.h"
#include "roots.h"
#include <math.h>
#include <stdlib.h>

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

// den_b with the factor (x - e) removed at each end e where U vanishes; false when that factor
// is squared there (b·y ~ (x - e)^(1/2 - k) is integrable only for k <= 1)
static bool strip_ends(const SqrtShape* s, const Rational* const* ends, Polynomial** rest) {
    *rest = poly_copy(s->den_b);
    bool ok = true;
    Rational v;
    rat_init(&v);
    for (int i = 0; i < 2 && ok && s->U; i++) {
        poly_eval_rat(&v, s->U, ends[i]);
        if (!rat_is_zero(&v)) continue;
        Rational one = rat_one_val(), m;
        rat_init(&m);
        rat_neg(&m, ends[i]);
        Polynomial* lin = create_polynomial();
        add_term(lin, &one, 1); add_term(lin, &m, 0);
        for (int k = 0; ok; k++) {
            PolyDivResult d = poly_divmod(*rest, lin);
            bool divides = poly_is_zero(d.remainder);
            free_polynomial(d.remainder);
            if (!divides) { free_polynomial(d.quotient); break; }
            ok = k == 0;
            free_polynomial(*rest);
            *rest = d.quotient;
        }
        free_polynomial(lin);
        rat_clear(&one); rat_clear(&m);
    }
    rat_clear(&v);
    return ok;
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

// arg's value at x0: through Q(x)(y) (with limits), else substituted and folded by expr_at, which
// multiplies constant roots into one (sqrt(3)·sqrt(U(x0))); NULL when neither gives a value
static Expr* value_of(const Expr* arg, const At* p) {
    if (arg->type == EXPR_BINOP && (arg->binop.op == OP_MUL || arg->binop.op == OP_DIV)
        && (!has_x(arg->binop.right) || (arg->binop.op == OP_MUL && !has_x(arg->binop.left)))) {
        bool right = !has_x(arg->binop.right);                // g·c, g/c or c·g: c x-free (sqrt(2)), g's limit
        Expr* g = value_of(right ? arg->binop.left : arg->binop.right, p);
        if (!g) return NULL;
        Expr* c = expr_copy(right ? arg->binop.right : arg->binop.left);
        Expr* r = arg->binop.op == OP_DIV ? create_div(g, c) : right ? create_mul(g, c) : create_mul(c, g);
        Expr* f = expr_at(r, p->x0);
        free_expr(r);
        return f;
    }
    Expr* v = sqrt_value_at(arg, p->U, p->x0);
    if (v) return v;
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

// F(b) - F(a), each term evaluated and folded; NULL when a term has no value at an end
static Expr* exact_difference(const Expr* F, const Polynomial* U, const Rational* const* ends) {
    Expr** T = NULL;
    int n = 0;
    summands(F, 1, &T, &n);
    Expr** parts = malloc((size_t)(2 * n + 1) * sizeof(Expr*));
    int m = 0;
    bool ok = !expr_has_rootsum(F);
    for (int i = 0; i < n && ok; i++)
        for (int e = 0; e < 2 && ok; e++) {                      // + F_i(b), - F_i(a)
            const Rational* x0 = ends[1 - e];
            Expr* v = term_at(T[i], &(At){ U, x0, x0 == ends[0] ? (rat_cmp(ends[0], ends[1]) < 0 ? 1 : -1) : (rat_cmp(ends[1], ends[0]) < 0 ? 1 : -1) });
            ok = v != NULL;
            if (ok) {
                Expr* f = expr_at(v, ends[1 - e]);
                free_expr(v);
                ok = f != NULL;
                if (ok) parts[m++] = e == 0 ? f : create_neg(f);
            }
        }
    for (int i = 0; i < n; i++) free_expr(T[i]);
    free(T);
    Expr* r = ok ? const_fold_sum(parts, m) : NULL;
    if (!ok) for (int i = 0; i < m; i++) free_expr(parts[i]);
    free(parts);
    return r;
}

// --- tanh-sinh quadrature, from each end ---

// f at e + s·u for u >= 0 into the range (s = +1 at the lower end, -1 at the upper): the
// polynomials of a + b·sqrt(U) shifted exactly to e, so u near 0 (a singular end) loses nothing
typedef struct {
    DPoly an, ad, bn, bd, U;
    bool root;
} Near;

// p(e + s·u) as a polynomial in u, exactly, in doubles
static DPoly shifted(const Polynomial* p, const Rational* e, int s) {
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
    DPoly d = dpoly_from(r);
    free_polynomial(r); free_polynomial(lin);
    rat_clear(&c); rat_clear(&sr);
    return d;
}

static Near near_end(const SqrtShape* s, const Rational* e, int side) {
    Near n = { shifted(s->num_a, e, side), shifted(s->den_a, e, side), shifted(s->num_b, e, side), shifted(s->den_b, e, side),
               { NULL, -1 }, s->U != NULL };
    if (s->U) n.U = shifted(s->U, e, side);
    return n;
}

static void near_free(Near* n) {
    dpoly_free(&n->an); dpoly_free(&n->ad); dpoly_free(&n->bn); dpoly_free(&n->bd);
    if (n->root) dpoly_free(&n->U);
}

static double near_at(const Near* n, double u) {
    double v = dpoly_eval(&n->an, u) / dpoly_eval(&n->ad, u);
    if (n->root) v += dpoly_eval(&n->bn, u) / dpoly_eval(&n->bd, u) * sqrt(fmax(dpoly_eval(&n->U, u), 0)); // U >= 0 here
    return v;
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

DefiniteResult integrate_definite_sqrt(const Expr* f, const Bounds* ab) {
    if (ab->a.inf || ab->b.inf || ab->a.pi || ab->b.pi) return fail("infinite bounds and multiples of pi with square roots are not supported yet");
    const Rational *a = ab->a.x, *b = ab->b.x;
    const Rational* lo = rat_cmp(a, b) <= 0 ? a : b, *hi = lo == a ? b : a;
    if (rat_cmp(a, b) == 0) return (DefiniteResult){ create_const(0), 0, NULL, NULL };
    SqrtShape s;
    const char* why = NULL;
    if (!sqrt_shape(f, &s, &why)) return fail(why);
    DefiniteResult r = { NULL, 0, NULL, NULL };
    const Rational* ends[2] = { lo, hi };
    Polynomial* rest = NULL;
    if (!real_on(s.U, lo, hi)) r.error = "the square root is not real on all of [a, b]";
    else if (root_in(s.den_a, &(Range){ lo, hi, true })) r.error = "the integrand has a pole in [a, b], so the integral diverges";
    else if (!strip_ends(&s, ends, &rest) || root_in(rest, &(Range){ lo, hi, true })) r.error = "the integrand has a pole in [a, b], so the integral diverges";
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
    free_polynomial(rest);
    sqrt_shape_free(&s);
    return r;
}
