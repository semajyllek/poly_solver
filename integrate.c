// Integration of polynomials and rational functions over Q.
//
//   ∫ N/D = ∫ q  +  g  +  ∫ a/D*        (polynomial part, Hermite reduction)
//
// Hermite reduction (Mack's linear version, Bronstein "Symbolic Integration I", §2.2)
// yields the rational part g without factoring; the remainder a/D* has a square-free
// denominator. That remainder is split by CRT over the Q-irreducible factors P of D*:
//   deg P = 1       -> c·ln|P|
//   deg P = 2       -> c·ln|P| plus atan (negative discriminant) or ln of a ratio (positive)
//   deg P >= 3      -> Lazard–Rioboo–Trager (logpart.c): logs/atans when the log coefficients
//                      have degree <= 2 over Q, otherwise an exact RootSum.
// Every result of the Hermite step is checked exactly (g' + a/D* == r/D) before use.
#include "expression.h"
#include "factor.h"
#include "logpart.h"
#include "rootform.h"
#include <stdlib.h>
#include <string.h>

static const char* last_error = NULL;

const char* integrate_error(void) { return last_error; }

typedef struct {
    Polynomial* num;
    Polynomial* den;
} RatFn;

// --- small polynomial helpers ---

static Polynomial* poly_one(void) {
    Rational one = rat_one_val();
    Polynomial* p = poly_constant(&one);
    rat_clear(&one);
    return p;
}

// exact quotient a / b (caller guarantees divisibility)
static Polynomial* divexact(const Polynomial* a, const Polynomial* b) {
    return poly_divide(a, b);
}

// (s, t) with s·a + t·b = c and deg s < deg b; requires gcd(a, b) = 1
typedef struct {
    Polynomial* s;
    Polynomial* t;
} Dioph;

static Dioph diophantine(const Polynomial* a, const Polynomial* b, const Polynomial* c) {
    PolyBezout bz = poly_xgcd(a, b);
    Polynomial* sc = poly_multiply(bz.s, c);
    Polynomial* s = poly_mod(sc, b);
    Polynomial* sa = poly_multiply(s, a);
    Polynomial* rest = poly_subtract(c, sa);
    Dioph d = { s, divexact(rest, b) };
    free_polynomial(sc); free_polynomial(sa); free_polynomial(rest);
    free_bezout(&bz);
    return d;
}

// acc += n/d, reduced by the gcd; consumes n and d
static void ratfn_add(RatFn* acc, Polynomial* n, Polynomial* d) {
    Polynomial* a = poly_multiply(acc->num, d);
    Polynomial* b = poly_multiply(n, acc->den);
    Polynomial* num = poly_add(a, b);
    Polynomial* den = poly_multiply(acc->den, d);
    Polynomial* g = poly_gcd(num, den);
    free_polynomial(acc->num); free_polynomial(acc->den);
    acc->num = poly_is_zero(num) ? create_polynomial() : divexact(num, g);
    acc->den = poly_is_zero(num) ? poly_one() : divexact(den, g);
    free_polynomial(a); free_polynomial(b); free_polynomial(num);
    free_polynomial(den); free_polynomial(g); free_polynomial(n); free_polynomial(d);
}

// --- Hermite reduction ---

typedef struct {
    RatFn g;          // rational part of the integral
    Polynomial* a;    // remaining integrand a / dstar
    Polynomial* dstar;
} Hermite;

// -Ds · Dm' / Dm, which is a polynomial
static Polynomial* hermite_coupling(const Polynomial* ds, const Polynomial* dm) {
    Polynomial* dmp = poly_derivative(dm);
    Polynomial* prod = poly_multiply(ds, dmp);
    Polynomial* q = divexact(prod, dm);
    Rational m1 = rat_from_int(-1);
    Polynomial* r = poly_scale(q, &m1);
    rat_clear(&m1);
    free_polynomial(dmp); free_polynomial(prod); free_polynomial(q);
    return r;
}

static Hermite hermite_reduce(const Polynomial* A0, const Polynomial* D) {
    Hermite h = { { create_polynomial(), poly_one() }, poly_copy(A0), NULL };
    Polynomial* dp = poly_derivative(D);
    Polynomial* dm = poly_gcd(D, dp);
    h.dstar = divexact(D, dm);
    free_polynomial(dp);
    while (polynomial_degree(dm) > 0) {
        Polynomial* dmp = poly_derivative(dm);
        Polynomial* dm2 = poly_gcd(dm, dmp);
        Polynomial* dms = divexact(dm, dm2);
        Polynomial* coupling = hermite_coupling(h.dstar, dm);
        Dioph bc = diophantine(coupling, dms, h.a);   // B·coupling + C·Dms = A
        Polynomial* bp = poly_derivative(bc.s);
        Polynomial* ratio = divexact(h.dstar, dms);
        Polynomial* corr = poly_multiply(bp, ratio);
        free_polynomial(h.a);
        h.a = poly_subtract(bc.t, corr);               // A <- C - B'·Ds/Dms
        ratfn_add(&h.g, bc.s, poly_copy(dm));          // g <- g + B/Dm
        free_polynomial(bc.t); free_polynomial(bp); free_polynomial(ratio); free_polynomial(corr);
        free_polynomial(dmp); free_polynomial(dms); free_polynomial(coupling); free_polynomial(dm);
        dm = dm2;
    }
    free_polynomial(dm);
    return h;
}

// exact check: g' + a/dstar == A0/D as rational functions
static bool hermite_ok(const Hermite* h, const Polynomial* A0, const Polynomial* D) {
    Polynomial *gn = h->g.num, *gd = h->g.den;
    Polynomial *gnp = poly_derivative(gn), *gdp = poly_derivative(gd);
    Polynomial *t1 = poly_multiply(gnp, gd), *t2 = poly_multiply(gn, gdp);
    Polynomial* dnum = poly_subtract(t1, t2);           // numerator of g' over gd^2
    Polynomial* gd2 = poly_multiply(gd, gd);
    Polynomial *u1 = poly_multiply(dnum, h->dstar), *u2 = poly_multiply(h->a, gd2);
    Polynomial* lnum = poly_add(u1, u2);                // lhs = lnum / (gd^2 · dstar)
    Polynomial* lden = poly_multiply(gd2, h->dstar);
    Polynomial *x = poly_multiply(lnum, D), *y = poly_multiply(A0, lden);
    Polynomial* diff = poly_subtract(x, y);
    bool ok = poly_is_zero(diff);
    Polynomial* parts[] = { gnp, gdp, t1, t2, dnum, gd2, u1, u2, lnum, lden, x, y, diff };
    for (size_t i = 0; i < sizeof parts / sizeof parts[0]; i++) free_polynomial(parts[i]);
    return ok;
}

static void hermite_free(Hermite* h) {
    free_polynomial(h->g.num); free_polynomial(h->g.den);
    free_polynomial(h->a); free_polynomial(h->dstar);
}

// --- expression builders ---

// acc + t, flattening t's own additions so signs print as "a - b" rather than "a + -b"
static Expr* sum(Expr* acc, Expr* t) {
    if (!t) return acc;
    if (!acc) return t;
    if (t->type == EXPR_BINOP && t->binop.op == OP_ADD) {
        Expr *l = t->binop.left, *r = t->binop.right;
        free(t); // only the shell; its children are re-linked below
        return sum(sum(acc, l), r);
    }
    return create_add(acc, t);
}

static Expr* scaled(const Rational* c, Expr* e) {
    return create_mul(create_const_rat(c), e);
}

static Expr* sqrt_of(const mpz_t m) {
    Rational r;
    rat_init(&r);
    mpq_set_z(r.val, m);
    Expr* e = create_sqrt(create_const_rat(&r));
    rat_clear(&r);
    return e;
}

// num/den with a primitive positive denominator; a leading minus becomes -1 · (...)
static Expr* ratfn_expr(const RatFn* g) {
    Rational c = poly_content(g->den), inv;
    rat_init(&inv);
    if (rat_is_negative(&g->den->head->coeff)) rat_neg(&c, &c);
    mpq_inv(inv.val, c.val);
    Polynomial* den = poly_scale(g->den, &inv);
    Polynomial* num = poly_scale(g->num, &inv);
    rat_clear(&c); rat_clear(&inv);
    if (polynomial_degree(den) == 0) { free_polynomial(den); return create_poly_expr(num); }
    bool neg = rat_is_negative(&num->head->coeff);
    if (neg) {
        Rational m1 = rat_from_int(-1);
        Polynomial* pos = poly_scale(num, &m1);
        free_polynomial(num);
        num = pos;
        rat_clear(&m1);
    }
    Expr* e = create_div(create_poly_expr(num), create_poly_expr(den));
    if (!neg) return e;
    Rational m1 = rat_from_int(-1);
    e = scaled(&m1, e);
    rat_clear(&m1);
    return e;
}

// --- integrals of num / P for one irreducible P (deg num < deg P) ---

static Expr* ln_term(const Rational* c, const Polynomial* P) {
    if (rat_is_zero(c)) return NULL;
    return scaled(c, create_ln(create_poly_expr(poly_copy(P))));
}

// (2a x + b) / r as a polynomial
static Polynomial* shifted_linear(const Rational* a, const Rational* b, const Rational* r) {
    Polynomial* p = create_polynomial();
    Rational t;
    rat_init(&t);
    rat_div(&t, a, r);
    mpq_mul_2exp(t.val, t.val, 1);
    add_term(p, &t, 1);
    rat_div(&t, b, r);
    add_term(p, &t, 0);
    rat_clear(&t);
    return p;
}

// K · ∫ dx / (a x^2 + b x + c) for an irreducible quadratic
static Expr* quad_reciprocal(const Rational* K, const Polynomial* P) {
    Rational a, b, c, d, t, r;
    mpz_t m;
    rat_init(&a); rat_init(&b); rat_init(&c); rat_init(&d); rat_init(&t); rat_init(&r);
    mpz_init(m);
    poly_coeff(&a, P, 2); poly_coeff(&b, P, 1); poly_coeff(&c, P, 0);
    rat_mul(&d, &b, &b);                      // D = b^2 - 4ac
    rat_mul(&t, &a, &c);
    mpq_mul_2exp(t.val, t.val, 2);
    rat_sub(&d, &d, &t);
    bool neg_disc = rat_is_negative(&d);
    rat_abs(&d, &d);
    rat_sqrt_split(&d, &r, m);                // sqrt|D| = r·sqrt(m)
    Polynomial* lin = shifted_linear(&a, &b, &r);
    bool unit_m = mpz_cmp_ui(m, 1) == 0;
    // coefficient: 2K/(r·sqrt m) for atan, K/(r·sqrt m) for the ln ratio; rationalized as /(r·m)·sqrt m
    rat_div(&t, K, &r);
    if (neg_disc) mpq_mul_2exp(t.val, t.val, 1);
    if (!unit_m) { mpq_set_z(d.val, m); rat_div(&t, &t, &d); }
    Expr* e;
    if (neg_disc) {
        Expr* arg = unit_m ? create_poly_expr(lin) : create_div(create_poly_expr(lin), sqrt_of(m));
        e = create_atan(arg);
    } else { // ln|(lin - sqrt m)/(lin + sqrt m)|
        Expr* lo = create_sub(create_poly_expr(poly_copy(lin)), sqrt_of(m));
        Expr* hi = create_add(create_poly_expr(lin), sqrt_of(m));
        e = create_ln(create_div(lo, hi));
    }
    if (!unit_m) e = create_mul(sqrt_of(m), e);
    e = scaled(&t, e);
    rat_clear(&a); rat_clear(&b); rat_clear(&c); rat_clear(&d); rat_clear(&t); rat_clear(&r);
    mpz_clear(m);
    return e;
}

// ∫ (B x + C) / P = B/(2a)·ln|P| + (C - B·b/(2a))·∫ 1/P
static Expr* integrate_quadratic(const Polynomial* num, const Polynomial* P) {
    Rational a, b, B, C, lc, K;
    rat_init(&a); rat_init(&b); rat_init(&B); rat_init(&C); rat_init(&lc); rat_init(&K);
    poly_coeff(&a, P, 2); poly_coeff(&b, P, 1);
    poly_coeff(&B, num, 1); poly_coeff(&C, num, 0);
    rat_div(&lc, &B, &a);
    mpq_div_2exp(lc.val, lc.val, 1);          // B / (2a)
    rat_mul(&K, &lc, &b);
    rat_sub(&K, &C, &K);                      // C - B·b/(2a)
    Expr* e = ln_term(&lc, P);
    if (!rat_is_zero(&K)) e = sum(e, quad_reciprocal(&K, P));
    rat_clear(&a); rat_clear(&b); rat_clear(&B); rat_clear(&C); rat_clear(&lc); rat_clear(&K);
    return e;
}

static bool expand_roots = true; // integrate_rootsum turns this off

static bool has_rootsum(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_ROOTSUM) return true;
    return e->type == EXPR_BINOP && (has_rootsum(e->binop.left) || has_rootsum(e->binop.right));
}

static Expr* integrate_simple(const Polynomial* num, const Polynomial* P) {
    if (poly_is_zero(num)) return NULL;
    int d = polynomial_degree(P);
    if (d == 1) {
        Rational k, a;
        rat_init(&k); rat_init(&a);
        poly_coeff(&k, num, 0);
        poly_coeff(&a, P, 1);
        rat_div(&k, &k, &a);                  // ∫ k/(a x + b) = k/a · ln|a x + b|
        Expr* e = ln_term(&k, P);
        rat_clear(&k); rat_clear(&a);
        return e;
    }
    if (d == 2) return integrate_quadratic(num, P);
    Expr* e = lrt_log_part(num, P);                      // Lazard–Rioboo–Trager
    if (!expand_roots || !has_rootsum(e)) return e;
    Expr* r = root_log_part(num, P);                     // RootSum -> radicals when P's roots are known
    if (!r) return e;
    free_expr(e);
    return r;
}

// ∫ a / dstar with dstar square-free: CRT split over the irreducible factors
static Expr* log_part(const Polynomial* a, const Polynomial* dstar) {
    Factorization* F = factorize(dstar);
    Expr* result = NULL;
    for (int i = 0; i < F->count && !last_error; i++) {
        Polynomial* q = poly_constant(&F->content); // cofactor content · Π_{j≠i} P_j
        for (int j = 0; j < F->count; j++) {
            if (j == i) continue;
            Polynomial* t = poly_multiply(q, F->factors[j]);
            free_polynomial(q);
            q = t;
        }
        Dioph st = diophantine(q, F->factors[i], a); // s·Q + t·P = a  =>  a/(QP) = s/P + t/Q
        result = sum(result, integrate_simple(st.s, F->factors[i]));
        free_polynomial(st.s); free_polynomial(st.t); free_polynomial(q);
    }
    free_factorization(F);
    if (last_error) { free_expr(result); return NULL; }
    return result;
}

// --- final cleanup: merge c1·f + c2·f into (c1+c2)·f ---

typedef struct {
    Expr** terms;
    int count;
} TermList;

static void collect_terms(Expr* e, TermList* L) { // consumes e
    if (e->type == EXPR_BINOP && e->binop.op == OP_ADD) {
        collect_terms(e->binop.left, L);
        collect_terms(e->binop.right, L);
        free(e);
        return;
    }
    L->terms = realloc(L->terms, (size_t)(L->count + 1) * sizeof(Expr*));
    L->terms[L->count++] = e;
}

// for c·f return f's printed form (the merge key), else NULL
static char* scaled_key(const Expr* e) {
    if (e->type != EXPR_BINOP || e->binop.op != OP_MUL || e->binop.left->type != EXPR_CONST) return NULL;
    return expr_to_string(e->binop.right);
}

// fold term j into term i when both are constant multiples of the same expression
static bool merge_into(Expr* ti, Expr* tj) {
    char *ki = scaled_key(ti), *kj = scaled_key(tj);
    bool same = ki && kj && strcmp(ki, kj) == 0;
    if (same) rat_add(&ti->binop.left->constant, &ti->binop.left->constant, &tj->binop.left->constant);
    free(ki); free(kj);
    return same;
}

static Expr* combine_like_terms(Expr* e) {
    TermList L = { NULL, 0 };
    collect_terms(e, &L);
    Expr* out = NULL;
    for (int i = 0; i < L.count; i++) {
        if (!L.terms[i]) continue;
        for (int j = i + 1; j < L.count; j++)
            if (L.terms[j] && merge_into(L.terms[i], L.terms[j])) { free_expr(L.terms[j]); L.terms[j] = NULL; }
        bool zero = L.terms[i]->type == EXPR_BINOP && L.terms[i]->binop.op == OP_MUL &&
                    L.terms[i]->binop.left->type == EXPR_CONST && rat_is_zero(&L.terms[i]->binop.left->constant);
        if (zero) free_expr(L.terms[i]);
        else out = sum(out, L.terms[i]);
    }
    free(L.terms);
    return out ? out : create_poly_expr(create_polynomial());
}

// --- entry points ---

static Expr* integrate_rational(const Polynomial* num, const Polynomial* den) {
    PolyDivResult dr = poly_divmod(num, den);
    Expr* result = poly_is_zero(dr.quotient) ? NULL : create_poly_expr(poly_integral(dr.quotient));
    if (!poly_is_zero(dr.remainder)) {
        Hermite h = hermite_reduce(dr.remainder, den);
        if (!hermite_ok(&h, dr.remainder, den)) last_error = "internal error: Hermite reduction check failed";
        if (!last_error && !poly_is_zero(h.g.num)) result = sum(result, ratfn_expr(&h.g));
        if (!last_error && polynomial_degree(h.a) >= polynomial_degree(h.dstar)) {
            PolyDivResult extra = poly_divmod(h.a, h.dstar); // keep the log part proper
            result = sum(result, create_poly_expr(poly_integral(extra.quotient)));
            free_polynomial(extra.quotient);
            free_polynomial(h.a);
            h.a = extra.remainder;
        }
        if (!last_error && !poly_is_zero(h.a)) result = sum(result, log_part(h.a, h.dstar));
        hermite_free(&h);
    }
    free_polynomial(dr.quotient);
    free_polynomial(dr.remainder);
    if (last_error) { free_expr(result); return NULL; }
    return result ? combine_like_terms(result) : create_poly_expr(create_polynomial());
}

Expr* integrate(Expr* expr, const char* variable) {
    (void)variable; // single variable for now
    last_error = NULL;
    if (!expr) { last_error = "no expression"; return NULL; }
    Expr* canonical = canonicalize(expr);
    Expr* result = NULL;
    if (!canonical) last_error = "only polynomials and rational functions can be integrated";
    else if (canonical->type == EXPR_POLY) result = create_poly_expr(poly_integral(canonical->poly));
    else if (canonical->type == EXPR_RATIONAL_FN)
        result = integrate_rational(canonical->ratfn.num, canonical->ratfn.den);
    free_expr(canonical);
    return result;
}

Expr* integrate_rootsum(Expr* expr, const char* variable) {
    expand_roots = false;
    Expr* r = integrate(expr, variable);
    expand_roots = true;
    return r;
}
