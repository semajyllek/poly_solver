// Real closed forms for ∫ a/P when the roots of the irreducible P are known exactly:
//   binomial    P = p_n·x^n + p_0          roots ρ·e^(iπ·j/n), ρ = |p_0/p_n|^(1/n)
//   cyclotomic  P = Q(x^s), Q = q_m·b^m·Φ_N(y/b), b ∈ Q
//               roots b^(1/s)·e^(iπ·(2k + 2lN)/(sN)), gcd(k, N) = 1, 0 <= l < s
// With h = a·(P')^(-1) mod P we have ∫ a/P = Σ_r h(r)·ln(x - r). A real root r gives
// h(r)·ln|x - r|; a conjugate pair r = ρ·e^(±iθ), 0 < θ < π, with h(r) = H_re + i·H_im gives
//   H_re·ln(x² - 2ρ·cosθ·x + ρ²) - 2·H_im·atan((x - ρ·cosθ) / (ρ·sinθ)),
// whose atan argument never has a pole, so no Rioboo recursion is needed.
#include "rootform.h"
#include "surd.h"
#include <stdlib.h>

// the roots are R^(1/n)·e^(iπ·j/L) for j in js
typedef struct {
    Rational R;
    int n;
    int L;
    int* js;
    int count;
} RootSet;

// ρ = R^(1/n) is set by the caller; room for `roots` angles
static void root_set_init(RootSet* rs, const Rational* R, int roots) {
    rat_init_set(&rs->R, R);
    rs->js = malloc((size_t)roots * sizeof(int));
    rs->count = 0;
}

static void root_set_free(RootSet* rs) {
    rat_clear(&rs->R);
    free(rs->js);
}

// --- recognizing P ---

// P = p_n·x^n + p_0: x^n = c with c = -p_0/p_n, roots |c|^(1/n)·e^(iπ(2k + [c < 0])/n)
static bool binomial_roots(const Polynomial* P, RootSet* rs) {
    const Term* lead = P->head;
    if (!lead || !lead->next || lead->next->next || lead->next->exponent != 0) return false;
    int n = lead->exponent;
    Rational c;
    rat_init(&c);
    rat_div(&c, &lead->next->coeff, &lead->coeff);
    rat_neg(&c, &c);
    bool neg = rat_is_negative(&c);
    rat_abs(&c, &c);
    root_set_init(rs, &c, n);
    rs->n = rs->L = n;
    for (int k = 0; k < n; k++) rs->js[rs->count++] = 2 * k + neg;
    rat_clear(&c);
    return true;
}

static int mobius(int n) {
    int r = 1;
    for (int p = 2; p * p <= n; p++)
        if (n % p == 0) {
            n /= p;
            if (n % p == 0) return 0;
            r = -r;
        }
    return n > 1 ? -r : r;
}

static int gcd_i(int a, int b) { return b ? gcd_i(b, a % b) : a; }

static int euler_phi(int n) {
    int r = 0;
    for (int k = 1; k <= n; k++) r += gcd_i(k, n) == 1;
    return r;
}

// Φ_N = Π_{d | N} (x^d - 1)^μ(N/d)
static Polynomial* cyclotomic(int N) {
    Rational one = rat_one_val(), m1 = rat_from_int(-1);
    Polynomial *num = poly_constant(&one), *den = poly_constant(&one);
    for (int d = 1; d <= N; d++) {
        int mu = (N % d == 0) ? mobius(N / d) : 0;
        if (mu == 0) continue;
        Polynomial* t = create_polynomial();
        add_term(t, &one, d);
        add_term(t, &m1, 0);
        Polynomial** side = mu > 0 ? &num : &den;
        Polynomial* prod = poly_multiply(*side, t);
        free_polynomial(*side);
        free_polynomial(t);
        *side = prod;
    }
    Polynomial* phi = poly_divide(num, den);
    free_polynomial(num); free_polynomial(den);
    rat_clear(&one); rat_clear(&m1);
    return phi;
}

// the positive rational b with b^n = q, if there is one
static bool rational_root(Rational* b, const Rational* q, int n) {
    if (rat_is_negative(q) || rat_is_zero(q)) return false;
    bool exact = mpz_root(mpq_numref(b->val), mpq_numref(q->val), (unsigned long)n)
              && mpz_root(mpq_denref(b->val), mpq_denref(q->val), (unsigned long)n);
    mpq_canonicalize(b->val);
    return exact;
}

// P == p_n·b^(n-k)·Φ_N[k] coefficient by coefficient
static bool matches_scaled(const Polynomial* P, const Polynomial* phi, const Rational* b) {
    int n = polynomial_degree(P);
    Rational pn, pk, fk, want, bk;
    rat_init(&pn); rat_init(&pk); rat_init(&fk); rat_init(&want); rat_init(&bk);
    poly_coeff(&pn, P, n);
    bool ok = polynomial_degree(phi) == n;
    for (int k = 0; ok && k <= n; k++) {
        poly_coeff(&pk, P, k);
        poly_coeff(&fk, phi, k);
        rat_pow_int(&bk, b, n - k);
        rat_mul(&want, &pn, &bk);
        rat_mul(&want, &want, &fk);
        ok = rat_eq(&want, &pk);
    }
    rat_clear(&pn); rat_clear(&pk); rat_clear(&fk); rat_clear(&want); rat_clear(&bk);
    return ok;
}

// Q = q_m·b^m·Φ_N(y/b) for some N with φ(N) = m and rational b > 0: N, or 0
static int scaled_cyclotomic(const Polynomial* Q, Rational* b) {
    int m = polynomial_degree(Q), found = 0;
    Rational q, q0, qm;
    rat_init(&q); rat_init(&q0); rat_init(&qm);
    poly_coeff(&q0, Q, 0);
    poly_coeff(&qm, Q, m);
    rat_div(&q, &q0, &qm);                                // b^m, since Φ_N(0) = 1 for N >= 2
    if (rational_root(b, &q, m))
        for (int N = 3; N <= 2 * m * m + 2 && !found; N++) { // φ(N) >= sqrt(N/2)
            if (euler_phi(N) != m) continue;
            Polynomial* phi = cyclotomic(N);
            if (matches_scaled(Q, phi, b)) found = N;
            free_polynomial(phi);
        }
    rat_clear(&q); rat_clear(&q0); rat_clear(&qm);
    return found;
}

// gcd of the exponents of P, so that P(x) = Q(x^g)
static int exponent_gcd(const Polynomial* P) {
    int g = 0;
    for (const Term* t = P->head; t; t = t->next) g = gcd_i(g, t->exponent);
    return g;
}

// Q with Q(x^s) = P
static Polynomial* deflate(const Polynomial* P, int s) {
    Polynomial* Q = create_polynomial();
    for (const Term* t = P->head; t; t = t->next) add_term(Q, &t->coeff, t->exponent / s);
    return Q;
}

// P(x) = Q(x^s) with Q a rescaled cyclotomic of degree >= 2: every s-th root of every b·ζ_N^k
static bool cyclotomic_roots(const Polynomial* P, RootSet* rs) {
    int g = exponent_gcd(P), N = 0;
    for (int s = 1; s <= g && !N; s++) {
        if (g % s) continue;
        Polynomial* Q = deflate(P, s);
        Rational b;
        rat_init(&b);
        N = polynomial_degree(Q) >= 2 ? scaled_cyclotomic(Q, &b) : 0;
        if (N) {
            root_set_init(rs, &b, polynomial_degree(P));
            rs->n = s;
            rs->L = s * N;
            for (int l = 0; l < s; l++)
                for (int k = 1; k < N; k++)
                    if (gcd_i(k, N) == 1) rs->js[rs->count++] = 2 * k + 2 * l * N;
        }
        rat_clear(&b);
        free_polynomial(Q);
    }
    return N != 0;
}

// --- residues and their real and imaginary parts ---

// h = a·(P')^(-1) mod P, so that h(r) = a(r)/P'(r) at every root r
static Polynomial* residues(const Polynomial* a, const Polynomial* P) {
    Polynomial* dp = poly_derivative(P);
    PolyBezout bz = poly_xgcd(dp, P);                     // s·P' + t·P = 1
    Polynomial* prod = poly_multiply(a, bz.s);
    Polynomial* h = poly_mod(prod, P);
    free_polynomial(prod); free_polynomial(dp);
    free_bezout(&bz);
    return h;
}

// one root ρ·e^(iπθ) of the set, with the residue polynomial h
typedef struct {
    const RootSet* rs;
    const Polynomial* h;
    Rational theta;
} RootAt;

// c·ρ^k·trig(kθ·mult)
typedef struct {
    long c;
    int k;
    TrigKind kind;
} RhoSpec;

// c·ρ^k·trig(θ) as a one-monomial constant (no trig for TRIG_NONE)
static RealConst rho_term(const RootAt* at, RhoSpec sp) {
    Rational c = rat_from_int(sp.c), e = rat_from_si(sp.k, at->rs->n);
    Mono m;
    mono_init_radical(&m, &c, &(Power){ &at->rs->R, &e });
    if (sp.kind != TRIG_NONE) mono_mul_trig(&m, sp.kind, &at->theta);
    RealConst k = rc_zero();
    rc_add(&k, &m);
    rat_clear(&c); rat_clear(&e);
    return k;
}

// scale·Re h(r) (cos) or scale·Im h(r) (sin): Σ_k scale·h_k·ρ^k·trig(kθ)
static RealConst h_part(const RootAt* at, TrigKind kind, long scale) {
    RealConst sum = rc_zero();
    Rational c, e, ang, s = rat_from_int(scale);
    rat_init(&c); rat_init(&e); rat_init(&ang);
    for (const Term* t = at->h->head; t; t = t->next) {
        rat_mul(&c, &t->coeff, &s);
        rat_set_si(&e, t->exponent, at->rs->n);
        Rational k = rat_from_int(t->exponent);
        rat_mul(&ang, &k, &at->theta);
        rat_clear(&k);
        Mono m;
        mono_init_radical(&m, &c, &(Power){ &at->rs->R, &e });
        mono_mul_trig(&m, kind, &ang);
        rc_add(&sum, &m);
    }
    rat_clear(&c); rat_clear(&e); rat_clear(&ang); rat_clear(&s);
    return sum;
}

static Expr* add(Expr* acc, Expr* t) {
    if (!acc) return t;
    return t ? create_add(acc, t) : acc;
}

// h(r)·ln|x - r| for r = ±ρ
static Expr* real_root_term(const RootAt* at) {
    RealConst H = h_part(at, TRIG_COS, 1), K = rho_term(at, (RhoSpec){ -1, 1, TRIG_COS });
    Expr* e = rc_scale(&H, create_ln(rc_add_terms(rc_x_pow(1), &K, -1)));
    rc_clear(&H); rc_clear(&K);
    return e;
}

// H_re·ln(x² - 2ρcosθ·x + ρ²) - 2·H_im·atan(cscθ/ρ·x - cotθ) for the pair ρ·e^(±iθ)
static Expr* pair_term(const RootAt* at) {
    RealConst Hre = h_part(at, TRIG_COS, 1), Him = h_part(at, TRIG_SIN, -2);
    RealConst K1 = rho_term(at, (RhoSpec){ -2, 1, TRIG_COS }), K0 = rho_term(at, (RhoSpec){ 1, 2, TRIG_NONE });
    RealConst K2 = rho_term(at, (RhoSpec){ 1, -1, TRIG_CSC }), K3 = rho_term(at, (RhoSpec){ -1, 0, TRIG_COT });
    Expr* ln_arg = rc_add_terms(rc_add_terms(rc_x_pow(2), &K1, 1), &K0, -1);
    Expr* atan_arg = rc_add_terms(rc_add_terms(NULL, &K2, 1), &K3, -1);
    Expr* e = add(rc_scale(&Hre, create_ln(ln_arg)), rc_scale(&Him, create_atan(atan_arg)));
    rc_clear(&Hre); rc_clear(&Him); rc_clear(&K1); rc_clear(&K0); rc_clear(&K2); rc_clear(&K3);
    return e;
}

Expr* root_log_part(const Polynomial* a, const Polynomial* P) {
    RootSet rs;
    if (polynomial_degree(P) < 3 || (!binomial_roots(P, &rs) && !cyclotomic_roots(P, &rs))) return NULL;
    Polynomial* h = residues(a, P);
    Rational one = rat_one_val();
    Expr* result = NULL;
    for (int i = 0; i < rs.count; i++) {
        RootAt at = { &rs, h, rat_from_si(rs.js[i], rs.L) }; // θ in [0, 2), units of π
        int c = rat_cmp(&at.theta, &one);
        if (rat_is_zero(&at.theta) || c == 0) result = add(result, real_root_term(&at));
        else if (c < 0) result = add(result, pair_term(&at)); // θ > 1 is the conjugate of 2 - θ
        rat_clear(&at.theta);
    }
    rat_clear(&one);
    free_polynomial(h);
    root_set_free(&rs);
    return result;
}
