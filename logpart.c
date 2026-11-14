// Lazard–Rioboo–Trager: the logarithmic part of ∫ a/d for square-free d, deg a < deg d.
// (Bronstein, "Symbolic Integration I", §2.5 IntRationalLogPart and §2.8 LogToReal.)
//
//   R(t) = res_x(d, a - t·d')  = Π Q_i(t)^i          (square-free decomposition)
//   S_i  = d if i = deg d, else the subresultant of degree i of (d, a - t·d'),
//          divided by gcd(A_j, Q_i)^j for the square-free parts A_j of lc_x(S_i)
//   ∫ a/d = Σ_i Σ_{Q_i(α)=0} α · log S_i(α, x)
//
// Each Q_i is factored over Q. A factor M of degree 1 gives α·ln|S(α,x)| with α rational;
// degree 2 gives a closed real form (logs with sqrt for real α; a log plus Rioboo's
// LogToAtan arctangents for complex α); degree >= 3 is returned as an exact RootSum.
#include "logpart.h"
#include "modgcd.h"
#include "subres.h"
#include "factor.h"
#include <stdlib.h>

// --- small helpers ---

static Expr* sum(Expr* acc, Expr* t) {
    if (!t) return acc;
    return acc ? create_add(acc, t) : t;
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

// primitive integer multiple of p with positive leading coefficient (same ln up to a constant)
static Polynomial* primitive(const Polynomial* p) {
    Polynomial* q = poly_primitive_part(p);
    if (q->head && rat_is_negative(&q->head->coeff)) {
        Rational m1 = rat_from_int(-1);
        Polynomial* r = poly_scale(q, &m1);
        free_polynomial(q);
        q = r;
        rat_clear(&m1);
    }
    return q;
}

// primitive form of a temporary, consuming it
static Polynomial* primitive_owned(Polynomial* p) {
    Polynomial* q = primitive(p);
    free_polynomial(p);
    return q;
}

// c·ln|P|, skipped when c = 0
static Expr* ln_term(const Rational* c, Polynomial* P) {
    if (rat_is_zero(c)) { free_polynomial(P); return NULL; }
    return scaled(c, create_ln(create_poly_expr(P)));
}

// --- S_i: choose and normalize ---

// s with s·a ≡ 1 (mod M): modular for large M (modgcd.c), else extended Euclid over Q
static Polynomial* inverse_mod(const Polynomial* a, const Polynomial* M) {
    Polynomial* s = polynomial_degree(M) >= 6 ? modular_inverse(a, M) : NULL;
    if (s) return s;
    PolyBezout bz = poly_xgcd(a, M);                     // s·a + t·M = 1
    s = poly_copy(bz.s);
    free_bezout(&bz);
    return s;
}

// divide every t-coefficient of S exactly by g
static void divide_coeffs(PolyTX* S, const Polynomial* g) {
    for (int k = 0; k <= S->deg; k++) {
        Polynomial* q = poly_divide(S->c[k], g);
        free_polynomial(S->c[k]);
        S->c[k] = q;
    }
}

// Lazard–Rioboo normalization: S <- S / Π gcd(A_j, Q)^j over the square-free parts of lc_x(S)
static void normalize_s(PolyTX* S, const Polynomial* Q) {
    SquareFreeDecomp lc = square_free_decompose(S->c[S->deg]);
    for (int j = 0; j < lc.count; j++) {
        Polynomial* g = poly_gcd_fast(lc.parts[j], Q);
        if (polynomial_degree(g) > 0) {
            Polynomial* gj = poly_pow(g, (unsigned)(j + 1));
            divide_coeffs(S, gj);
            free_polynomial(gj);
        }
        free_polynomial(g);
    }
    free_square_free(&lc);
}

// S reduced modulo M(t) and made monic in x (lc is invertible mod M after normalization)
static PolyTX reduce_monic(const PolyTX* S, const Polynomial* M) {
    PolyTX r = ptx_alloc(S->deg);
    for (int k = 0; k <= S->deg; k++) { free_polynomial(r.c[k]); r.c[k] = poly_mod(S->c[k], M); }
    while (r.deg > 0 && poly_is_zero(r.c[r.deg])) { free_polynomial(r.c[r.deg]); r.deg--; }
    Polynomial* s = inverse_mod(r.c[r.deg], M);        // s·lc ≡ 1 (mod M)
    for (int k = 0; k <= r.deg; k++) {
        Polynomial* prod = poly_multiply(r.c[k], s);
        free_polynomial(r.c[k]);
        r.c[k] = poly_mod(prod, M);
        free_polynomial(prod);
    }
    free_polynomial(s);
    return r;
}

// the x-polynomial formed by coefficient `power` of t in each S coefficient
static Polynomial* t_slice(const PolyTX* S, int power) {
    Polynomial* p = create_polynomial();
    Rational c;
    rat_init(&c);
    for (int k = 0; k <= S->deg; k++) { poly_coeff(&c, S->c[k], power); add_term(p, &c, k); }
    rat_clear(&c);
    return p;
}

// u + c·v
static Polynomial* add_scaled(const Polynomial* u, const Rational* c, const Polynomial* v) {
    Polynomial* cv = poly_scale(v, c);
    Polynomial* r = poly_add(u, cv);
    free_polynomial(cv);
    return r;
}

// --- degree 1: one rational alpha ---

static Expr* rational_root_term(const PolyTX* S, const Polynomial* M) {
    Rational alpha, m1;
    rat_init(&alpha); rat_init(&m1);
    poly_coeff(&alpha, M, 0);
    poly_coeff(&m1, M, 1);
    rat_div(&alpha, &alpha, &m1);
    rat_neg(&alpha, &alpha);                 // alpha = -M0/M1
    Polynomial* s = create_polynomial();
    Rational c;
    rat_init(&c);
    for (int k = 0; k <= S->deg; k++) { poly_eval_rat(&c, S->c[k], &alpha); add_term(s, &c, k); }
    Expr* e = ln_term(&alpha, primitive(s));
    free_polynomial(s);
    rat_clear(&c); rat_clear(&alpha); rat_clear(&m1);
    return e;
}

// --- degree 2: alpha = p ± r·sqrt(m) or p ± r·sqrt(m)·i, with S = S0 + t·S1 ---

typedef struct {
    Polynomial* U;  // S0 + p·S1
    Polynomial* W;  // r·S1
    QuadRoots q;
} QuadSplit;

static QuadSplit quad_split(const PolyTX* S, const Polynomial* M) {
    QuadSplit z;
    z.q = poly_quad_roots(M);
    Polynomial *s0 = t_slice(S, 0), *s1 = t_slice(S, 1);
    z.U = add_scaled(s0, &z.q.p, s1);
    z.W = poly_scale(s1, &z.q.r);
    free_polynomial(s0); free_polynomial(s1);
    return z;
}

static void quad_split_free(QuadSplit* z) {
    free_polynomial(z->U); free_polynomial(z->W); quad_roots_clear(&z->q);
}

// U^2 + sign·m·W^2
static Polynomial* norm_poly(const QuadSplit* z, int sign) {
    Polynomial *u2 = poly_multiply(z->U, z->U), *w2 = poly_multiply(z->W, z->W);
    Rational k;
    rat_init(&k);
    mpq_set_z(k.val, z->q.m);
    if (sign < 0) rat_neg(&k, &k);
    Polynomial* r = add_scaled(u2, &k, w2);
    free_polynomial(u2); free_polynomial(w2);
    rat_clear(&k);
    return r;
}

// sqrt(m)·W, printed as just sqrt(m) when W = 1
static Expr* sqrt_times_w(const QuadSplit* z) {
    bool unit = polynomial_degree(z->W) == 0 && rat_is_one(&z->W->head->coeff);
    return unit ? sqrt_of(z->q.m) : create_mul(sqrt_of(z->q.m), create_poly_expr(poly_copy(z->W)));
}

// U ± sqrt(m)·W as an expression
static Expr* u_pm_sqrt_w(const QuadSplit* z, bool plus) {
    Expr* w = sqrt_times_w(z);
    Expr* u = create_poly_expr(poly_copy(z->U));
    return plus ? create_add(u, w) : create_sub(u, w);
}

// real alpha: p·ln|U^2 - m·W^2| + r·sqrt(m)·ln|(U + sqrt(m)·W)/(U - sqrt(m)·W)|
static Expr* real_pair_terms(const QuadSplit* z) {
    Expr* e = ln_term(&z->q.p, primitive_owned(norm_poly(z, -1)));
    QuadSplit w = *z;                         // flip W to a positive leading coefficient:
    Rational c;                               // ln|(U+√m W)/(U-√m W)| is odd in W
    rat_init_set(&c, &z->q.r);
    bool flip = rat_is_negative(&z->W->head->coeff);
    if (flip) {
        Rational m1 = rat_from_int(-1);
        w.W = poly_scale(z->W, &m1);
        rat_neg(&c, &c);
        rat_clear(&m1);
    }
    Expr* ratio = create_ln(create_div(u_pm_sqrt_w(&w, true), u_pm_sqrt_w(&w, false)));
    e = sum(e, scaled(&c, create_mul(sqrt_of(z->q.m), ratio)));
    if (flip) free_polynomial(w.W);
    rat_clear(&c);
    return e;
}

// Rioboo's LogToAtan: Σ 2·atan(...) with derivative equal to that of i·log((A+iB)/(A-iB))
static Expr* log_to_atan(Polynomial* A, Polynomial* B, const Rational* coef) {
    Rational two;
    rat_init(&two);
    rat_add(&two, coef, coef);               // every term carries 2·coef
    Expr* e = NULL;
    while (!poly_is_zero(B)) {
        PolyDivResult dr = poly_divmod(A, B);
        bool divides = poly_is_zero(dr.remainder);
        free_polynomial(dr.remainder);
        if (divides) { e = sum(e, scaled(&two, create_atan(create_poly_expr(dr.quotient)))); break; }
        free_polynomial(dr.quotient);
        if (polynomial_degree(A) < polynomial_degree(B)) { // LogToAtan(A,B) = LogToAtan(-B, A)
            Rational m1 = rat_from_int(-1);
            Polynomial* nb = poly_scale(B, &m1);
            free_polynomial(B);
            B = A; A = nb;
            rat_clear(&m1);
            continue;
        }
        Rational m1 = rat_from_int(-1);
        Polynomial* negA = poly_scale(A, &m1);
        PolyBezout bz = poly_xgcd(B, negA);   // B·D - A·C = G
        Polynomial *ad = poly_multiply(A, bz.s), *bc = poly_multiply(B, bz.t);
        Polynomial* num = poly_add(ad, bc);
        Expr* arg = polynomial_degree(bz.g) == 0 ? create_poly_expr(poly_divide(num, bz.g))
                    : create_div(create_poly_expr(poly_copy(num)), create_poly_expr(poly_copy(bz.g)));
        e = sum(e, scaled(&two, create_atan(arg)));
        free_polynomial(A); free_polynomial(B);
        A = poly_copy(bz.s);
        B = poly_copy(bz.t);
        free_polynomial(ad); free_polynomial(bc); free_polynomial(num); free_polynomial(negA);
        free_bezout(&bz);
        rat_clear(&m1);
    }
    free_polynomial(A); free_polynomial(B);
    rat_clear(&two);
    return e;
}

// complex alpha = p ± r·sqrt(m)·i:
//   p·ln|U^2 + m·W^2| + r·LogToAtan(U, W)            when m = 1 (rational imaginary part)
//   p·ln|U^2 + m·W^2| - 2r·sqrt(m)·atan(sqrt(m)·W/U)  otherwise, or the equivalent
//   (up to a constant) 2r·sqrt(m)·atan(U/(sqrt(m)·W)) when U has the higher degree
static Expr* complex_pair_terms(const QuadSplit* z) {
    Expr* e = ln_term(&z->q.p, primitive_owned(norm_poly(z, +1)));
    if (mpz_cmp_ui(z->q.m, 1) == 0) // imaginary part r is rational; V = W = r·S1
        return sum(e, log_to_atan(poly_copy(z->U), poly_copy(z->W), &z->q.r));
    if (poly_is_zero(z->U)) return e;         // arg(i·W) is constant
    bool flip = polynomial_degree(z->U) > polynomial_degree(z->W); // arg = pi/2 - atan(U/V)
    Rational c;                               // ∓2·r (the sqrt(m) is a separate factor)
    rat_init(&c);
    rat_add(&c, &z->q.r, &z->q.r);
    if (!flip) rat_neg(&c, &c);
    Expr* arg = flip ? create_div(create_poly_expr(poly_copy(z->U)), sqrt_times_w(z))
                     : create_div(sqrt_times_w(z), create_poly_expr(poly_copy(z->U)));
    e = sum(e, scaled(&c, create_mul(sqrt_of(z->q.m), create_atan(arg))));
    rat_clear(&c);
    return e;
}

static Expr* quadratic_root_terms(const PolyTX* S, const Polynomial* M) {
    QuadSplit z = quad_split(S, M);
    Expr* e = z.q.is_complex ? complex_pair_terms(&z) : real_pair_terms(&z);
    quad_split_free(&z);
    return e;
}

// --- degree >= 3: exact RootSum ---

static Expr* rootsum_term(const PolyTX* S, const Polynomial* M) {
    Polynomial** s = malloc((size_t)(S->deg + 1) * sizeof(Polynomial*));
    for (int k = 0; k <= S->deg; k++) s[k] = poly_copy(S->c[k]);
    return create_rootsum(primitive(M), s, S->deg);
}

// all terms for the roots of one irreducible factor M of Q_i
static Expr* terms_for_factor(const PolyTX* S, const Polynomial* M) {
    PolyTX r = reduce_monic(S, M);
    int d = polynomial_degree(M);
    Expr* e = d == 1 ? rational_root_term(&r, M) : d == 2 ? quadratic_root_terms(&r, M) : rootsum_term(&r, M);
    ptx_free(&r);
    return e;
}

Expr* lrt_log_part(const Polynomial* a, const Polynomial* d) {
    RTPair pr = { d, a };
    SubresChain ch = rt_chain(&pr);
    SquareFreeDecomp sq = square_free_decompose(ch.S[0].c[0]); // the resultant R(t)
    Expr* result = NULL;
    for (int i = 1; i <= sq.count; i++) {
        const Polynomial* Q = sq.parts[i - 1];
        if (polynomial_degree(Q) <= 0) continue;
        PolyTX S = ptx_copy(&ch.S[i]); // d, g = a - t·d', or the subresultant S_i
        normalize_s(&S, Q);
        Factorization* F = factorize(Q);
        for (int f = 0; f < F->count; f++) result = sum(result, terms_for_factor(&S, F->factors[f]));
        free_factorization(F);
        ptx_free(&S);
    }
    free_square_free(&sq);
    chain_free(&ch);
    return result;
}
