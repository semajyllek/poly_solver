// Modular gcd and inverse over Q. Euclid over Q suffers coefficient blow-up (the remainders'
// numerators and denominators grow far beyond the answer's), which made reducing LRT's S_i
// modulo a factor of the resultant the slowest step for dense denominators. Here every
// computation is in F_p with word-size p; the answer's rational coefficients come back by
// CRT and rational reconstruction once enough primes agree, and are then checked exactly.
#include "modgcd.h"
#include "poly_zx.h"
#include <stdlib.h>

#define MAX_PRIMES 4000      // ~124,000 bits of modulus before giving up
#define FIRST_PRIME 2147483647UL // 2^31 - 1; zp arithmetic needs p < 2^31

// --- one modular problem: images, CRT, reconstruction, check ---

typedef struct {
    // the image for prime p into *out; returns its length (number of coefficients), or -1
    // when p is unusable (it divides a leading coefficient, or the image degenerates)
    int (*image)(const void* ctx, const ZpField* F, ZpPoly* out);
    // the candidate is the exact answer
    bool (*check)(const void* ctx, const Polynomial* candidate);
    const void* ctx;
} ModProblem;

// residues c[i] modulo m of the answer's coefficients
typedef struct {
    mpz_t* c;
    int len;
    mpz_t m;
} Crt;

static void crt_free(Crt* A) {
    for (int i = 0; i < A->len; i++) mpz_clear(A->c[i]);
    free(A->c);
    A->c = NULL;
    A->len = 0;
}

static unsigned long coeff_at(const ZpPoly* v, int i) { return i <= v->deg ? v->c[i] : 0; }

// start over from a single image (the first prime, or one with a smaller degree)
static void crt_reset(Crt* A, const ZpPoly* v, int len, unsigned long p) {
    crt_free(A);
    A->c = malloc((size_t)len * sizeof(mpz_t));
    A->len = len;
    for (int i = 0; i < len; i++) mpz_init_set_ui(A->c[i], coeff_at(v, i));
    mpz_set_ui(A->m, p);
}

// x ≡ c (mod m), x ≡ v (mod p): x = c + m·((v - c)·m⁻¹ mod p)
static void crt_add(Crt* A, const ZpPoly* v, unsigned long p) {
    mpz_t minv, t, pp;
    mpz_inits(minv, t, pp, NULL);
    mpz_set_ui(pp, p);
    mpz_invert(minv, A->m, pp);
    for (int i = 0; i < A->len; i++) {
        mpz_set_ui(t, coeff_at(v, i));
        mpz_sub(t, t, A->c[i]);
        mpz_mul(t, t, minv);
        mpz_mod(t, t, pp);
        mpz_addmul(A->c[i], A->m, t);
    }
    mpz_mul(A->m, A->m, pp);
    mpz_clears(minv, t, pp, NULL);
}

// n/d ≡ r (mod m) with |n|, d <= sqrt(m/2) (half-extended Euclid); false if there is none
static bool rat_reconstruct(Rational* out, const mpz_t r, const mpz_t m) {
    mpz_t r0, r1, t0, t1, q, tmp, bound;
    mpz_inits(r0, r1, t0, t1, q, tmp, bound, NULL);
    mpz_set(r0, m);
    mpz_mod(r1, r, m);
    mpz_set_ui(t0, 0);
    mpz_set_ui(t1, 1);
    mpz_fdiv_q_2exp(bound, m, 1);
    mpz_sqrt(bound, bound);
    while (mpz_cmp(r1, bound) > 0) {
        mpz_fdiv_q(q, r0, r1);
        mpz_submul(r0, q, r1); mpz_swap(r0, r1);          // (r0, r1) <- (r1, r0 - q·r1)
        mpz_submul(t0, q, t1); mpz_swap(t0, t1);          // (t0, t1) <- (t1, t0 - q·t1)
    }
    mpz_gcd(tmp, r1, t1);
    bool ok = mpz_sgn(t1) != 0 && mpz_cmpabs(t1, bound) <= 0 && mpz_cmp_ui(tmp, 1) == 0;
    if (ok) {
        mpq_set_num(out->val, r1);
        mpq_set_den(out->val, t1);                        // canonicalize fixes a negative d
        mpq_canonicalize(out->val);
    }
    mpz_clears(r0, r1, t0, t1, q, tmp, bound, NULL);
    return ok;
}

// the polynomial whose coefficients reconstruct from A, or NULL
static Polynomial* reconstruct(const Crt* A) {
    Polynomial* p = create_polynomial();
    Rational c;
    rat_init(&c);
    bool ok = true;
    for (int i = 0; ok && i < A->len; i++) {
        ok = rat_reconstruct(&c, A->c[i], A->m);
        if (ok) add_term(p, &c, i);
    }
    rat_clear(&c);
    if (!ok) { free_polynomial(p); return NULL; }
    return p;
}

static bool poly_same(const Polynomial* a, const Polynomial* b) {
    Polynomial* d = poly_subtract(a, b);
    bool same = poly_is_zero(d);
    free_polynomial(d);
    return same;
}

static unsigned long prev_prime(unsigned long p) {
    mpz_t z;
    mpz_init_set_ui(z, p);
    do mpz_sub_ui(z, z, 2); while (!mpz_probab_prime_p(z, 25));
    unsigned long q = mpz_get_ui(z);
    mpz_clear(z);
    return q;
}

// images for successive primes, keeping only those of the least length (unlucky primes give
// a longer gcd). Reconstruction is tried each time the prime count has grown by a quarter
// (so its total cost stays proportional to the last try), and a reconstruction that
// repeats at the next try is checked exactly.
static Polynomial* modular_solve(const ModProblem* pb) {
    Crt A = { .c = NULL, .len = 0 };
    mpz_init(A.m);
    Polynomial *last = NULL, *result = NULL;
    unsigned long p = FIRST_PRIME;
    int next_try = 1;
    for (int used = 0, tried = 0; !result && used < MAX_PRIMES && tried < 2 * MAX_PRIMES; p = prev_prime(p), tried++) {
        ZpField F = { p };
        ZpPoly v;
        int len = pb->image(pb->ctx, &F, &v);
        if (len < 0) continue;
        bool fresh = A.len == 0 || len < A.len;
        if (!fresh && len > A.len) { zp_free(&v); continue; }
        if (fresh) { crt_reset(&A, &v, len, p); free_polynomial(last); last = NULL; used = 0; next_try = 1; }
        else crt_add(&A, &v, p);
        zp_free(&v);
        if (++used < next_try) continue;
        next_try = used + (used / 4 > 1 ? used / 4 : 1);
        Polynomial* cand = reconstruct(&A);
        bool stable = cand && last && poly_same(cand, last);
        if (stable && pb->check(pb->ctx, cand)) result = cand;
        else { free_polynomial(last); last = cand; }
    }
    free_polynomial(last);
    crt_free(&A);
    mpz_clear(A.m);
    return result;
}

// --- integer forms ---

// c·p with integer coefficients, c the lcm of p's denominators
static Polynomial* integral_multiple(const Polynomial* p, Rational* c) {
    mpz_t l;
    mpz_init_set_ui(l, 1);
    for (const Term* t = p->head; t; t = t->next) mpz_lcm(l, l, mpq_denref(t->coeff.val));
    mpq_set_z(c->val, l);
    mpz_clear(l);
    return poly_scale(p, c);
}

// --- gcd ---

typedef struct {
    Polynomial *a, *b; // integer multiples of the inputs
    ZxPoly za, zb;
} GcdCtx;

static int gcd_image(const void* v, const ZpField* F, ZpPoly* out) {
    const GcdCtx* g = v;
    ZpPoly a = zx_to_zp(&g->za, F), b = zx_to_zp(&g->zb, F);
    int len = -1;
    if (a.deg == g->za.deg && b.deg == g->zb.deg) {      // p divides neither leading coefficient
        *out = zp_gcd(&a, &b, F);
        len = out->deg + 1;
    }
    zp_free(&a); zp_free(&b);
    return len;
}

static bool divides(const Polynomial* g, const Polynomial* a) {
    Polynomial* r = poly_mod(a, g);
    bool d = poly_is_zero(r);
    free_polynomial(r);
    return d;
}

// a monic common divisor has degree <= deg gcd, and no image has smaller degree than the
// gcd's, so a candidate of the least image degree dividing both is the gcd
static bool gcd_check(const void* v, const Polynomial* cand) {
    const GcdCtx* g = v;
    return polynomial_degree(cand) >= 0 && divides(cand, g->a) && divides(cand, g->b);
}

Polynomial* modular_gcd(const Polynomial* a, const Polynomial* b) {
    Rational ca, cb;
    rat_init(&ca); rat_init(&cb);
    GcdCtx g = { integral_multiple(a, &ca), integral_multiple(b, &cb), { NULL, -1 }, { NULL, -1 } };
    g.za = zx_from_poly(g.a);
    g.zb = zx_from_poly(g.b);
    ModProblem pb = { gcd_image, gcd_check, &g };
    Polynomial* r = modular_solve(&pb);
    zx_free(&g.za); zx_free(&g.zb);
    free_polynomial(g.a); free_polynomial(g.b);
    rat_clear(&ca); rat_clear(&cb);
    return r;
}

#define MODULAR_DEGREE 6

Polynomial* poly_gcd_fast(const Polynomial* a, const Polynomial* b) {
    bool big = polynomial_degree(a) >= MODULAR_DEGREE && polynomial_degree(b) >= MODULAR_DEGREE;
    Polynomial* g = big ? modular_gcd(a, b) : NULL;
    return g ? g : poly_gcd(a, b);
}

// --- inverse modulo M ---

typedef struct {
    Polynomial *a, *M; // integer multiples of the inputs (the same ideal for M)
    ZxPoly za, zM;
} InvCtx;

static int inv_image(const void* v, const ZpField* F, ZpPoly* out) {
    const InvCtx* c = v;
    ZpPoly a = zx_to_zp(&c->za, F), M = zx_to_zp(&c->zM, F);
    int len = -1;
    if (M.deg == c->zM.deg) {                             // p does not divide lc(M)
        ZpPoly g = zp_gcd(&a, &M, F);
        if (g.deg == 0) {                                 // a stays invertible modulo p
            ZpBezout bz = zp_xgcd(&a, &M, F);                 // deg s < deg M already
            *out = bz.s;
            zp_free(&bz.t);
            len = c->zM.deg;
        }
        zp_free(&g);
    }
    zp_free(&a); zp_free(&M);
    return len;
}

static bool inv_check(const void* v, const Polynomial* cand) {
    const InvCtx* c = v;
    Polynomial* prod = poly_multiply(cand, c->a);
    Polynomial* r = poly_mod(prod, c->M);
    bool one = polynomial_degree(r) == 0;
    if (one) {
        Rational k;
        rat_init(&k);
        poly_coeff(&k, r, 0);
        one = rat_is_one(&k);
        rat_clear(&k);
    }
    free_polynomial(prod); free_polynomial(r);
    return one;
}

Polynomial* modular_inverse(const Polynomial* a, const Polynomial* M) {
    Rational ca, cM;
    rat_init(&ca); rat_init(&cM);
    InvCtx c = { integral_multiple(a, &ca), integral_multiple(M, &cM), { NULL, -1 }, { NULL, -1 } };
    c.za = zx_from_poly(c.a);
    c.zM = zx_from_poly(c.M);
    ModProblem pb = { inv_image, inv_check, &c };
    Polynomial* s = modular_solve(&pb);                   // s·(ca·a) ≡ 1, so (ca·s)·a ≡ 1
    Polynomial* r = s ? poly_scale(s, &ca) : NULL;
    free_polynomial(s);
    zx_free(&c.za); zx_free(&c.zM);
    free_polynomial(c.a); free_polynomial(c.M);
    rat_clear(&ca); rat_clear(&cM);
    return r;
}
