// K1[θ2] and K1(θ2) for K1 = Q(x)(θ1): coefficient arithmetic in KFrac, the algorithms of
// kpoly_tmpl.h, and reduced fractions.
#include "k2poly.h"
#include <stdlib.h>

KFrac kf_zero(void) { return (KFrac){ { 0, 0, NULL }, kp_one() }; }
bool kf_is_zero(const KFrac* k) { return k->n.n == 0; }
KFrac kf_copy(const KFrac* k) { return (KFrac){ kp_copy(&k->n), kp_copy(&k->d) }; }

void k2_free(K2Poly* p) {
    for (int i = 0; i < p->n; i++) kfrac_free(&p->a[i]);
    free(p->a);
    p->a = NULL; p->n = 0;
}

K2Poly k2_zeros(int n, TowerKind kind) {
    K2Poly p = { malloc((size_t)(n > 0 ? n : 1) * sizeof(KFrac)), n, kind };
    for (int i = 0; i < n; i++) p.a[i] = kf_zero();
    return p;
}

void k2_trim(K2Poly* p) {
    while (p->n > 0 && kf_is_zero(&p->a[p->n - 1])) kfrac_free(&p->a[--p->n]);
}

K2Poly k2_const(KFrac c, TowerKind kind) {         // takes c
    K2Poly p = k2_zeros(1, kind);
    kfrac_free(&p.a[0]);
    p.a[0] = c;
    k2_trim(&p);
    return p;
}

// a[k] += c (taken)
void k2_accumulate(K2Poly* p, int k, KFrac c) {
    if (k >= p->n) {
        p->a = realloc(p->a, (size_t)(k + 1) * sizeof(KFrac));
        for (int i = p->n; i <= k; i++) p->a[i] = kf_zero();
        p->n = k + 1;
    }
    KFrac s = kfrac_add(&p->a[k], &c, p->kind);
    kfrac_free(&p->a[k]); kfrac_free(&c);
    p->a[k] = s;
}

// c·k
KFrac kf_scaled(const KFrac* k, const Rational* c, TowerKind kind) {
    (void)kind;
    if (rat_is_zero(c)) return kf_zero();
    QFn cq = qfn_rat(c);
    KFrac r = { kp_scale(&k->n, &cq), kp_copy(&k->d) };  // still reduced: c ≠ 0 is a unit
    qfn_free(&cq);
    return r;
}

K2Poly k2_add(const K2Poly* x, const K2Poly* y, int sign) {
    K2Poly r = k2_zeros(0, x->kind);
    Rational s = rat_from_int(sign);
    for (int i = 0; i < x->n; i++) k2_accumulate(&r, i, kf_copy(&x->a[i]));
    for (int i = 0; i < y->n; i++) k2_accumulate(&r, i, kf_scaled(&y->a[i], &s, x->kind));
    rat_clear(&s);
    k2_trim(&r);
    return r;
}

K2Poly k2_mul(const K2Poly* x, const K2Poly* y) {
    K2Poly r = k2_zeros(0, x->kind);
    for (int i = 0; i < x->n; i++)
        for (int j = 0; j < y->n; j++)
            if (!kf_is_zero(&x->a[i]) && !kf_is_zero(&y->a[j])) k2_accumulate(&r, i + j, kfrac_mul(&x->a[i], &y->a[j], x->kind));
    k2_trim(&r);
    return r;
}

int k2_deg(const K2Poly* p) { return p->n - 1; }

// --- the algorithms (kpoly_tmpl.h) over KFrac ---

// k·a
static KFrac kf_times_int(const KFrac* a, int k, TowerKind kind) {
    Rational kk = rat_from_int(k);
    KFrac r = kf_scaled(a, &kk, kind);
    rat_clear(&kk);
    return r;
}

// c·θ2^k (c copied)
static K2Poly k2_term(const KFrac* c, int k, TowerKind kind) {
    K2Poly p = k2_zeros(0, kind);
    if (!kf_is_zero(c)) k2_accumulate(&p, k, kf_copy(c));
    return p;
}

#define P K2Poly
#define C KFrac
#define CTX TowerKind
#define PFX(name) kt_##name
#define KDIV K2Div
#define KBEZ K2Bezout
#define KSF K2SquareFree
#define P_DEG(p) k2_deg(p)
#define P_LC(p) (&(p)->a[(p)->n - 1])
#define P_TERM(x, c, k) k2_term((c), (k), (x))
#define P_ADD(x, a, b, sign) k2_add((a), (b), (sign))
#define P_MUL(x, a, b) k2_mul((a), (b))
#define P_ONE(x) k2_const(kfrac_make(kp_one(), kp_one(), (x)), (x))
#define P_ZERO(x) k2_zeros(0, (x))
#define P_FREE(p) k2_free(p)
#define P_COEF(p, k) (&(p)->a[k])
#define P_LO(p) 0
#define P_N(p) ((p)->n)
#define C_MUL(x, a, b) kfrac_mul((a), (b), (x))
#define C_INV(x, a) kfrac_inv((a), (x))
#define C_FREE(a) kfrac_free(a)
#define C_SCALE_INT(x, a, k) kf_times_int((a), (k), (x))
#include "kpoly_tmpl.h"

K2Poly k2_copy(const K2Poly* a) { return kt_copy(a->kind, a); }
K2Poly k2_monic(const K2Poly* a) { return kt_monic(a->kind, a); }
K2Poly k2_dtheta(const K2Poly* a) { return kt_dtheta(a->kind, a); }
K2Div k2_divmod(const K2Poly* a, const K2Poly* b) { return kt_divmod(a->kind, a, b); }
K2Poly k2_exact_div(const K2Poly* a, const K2Poly* b) { return kt_exact_div(a->kind, a, b); }
K2Poly k2_gcd(const K2Poly* a, const K2Poly* b) { return kt_gcd(a->kind, a, b); }
K2Bezout k2_bezout(const K2Poly* a, const K2Poly* b, const K2Poly* c) { return kt_bezout(a->kind, a, b, c); }
K2SquareFree k2_square_free(const K2Poly* d) { return kt_square_free(d->kind, d); }

void k2div_free(K2Div* d) { k2_free(&d->q); k2_free(&d->r); }
void k2bezout_free(K2Bezout* z) { k2_free(&z->s); k2_free(&z->t); }

void k2sf_free(K2SquareFree* s) {
    for (int i = 0; i < s->m; i++) k2_free(&s->f[i]);
    free(s->f);
    s->f = NULL;
    s->m = 0;
}

// --- fractions ---

K2Frac k2frac_make(K2Poly n, K2Poly d) {
    K2Poly g = k2_gcd(&n, &d), nq = k2_exact_div(&n, &g), dq = k2_exact_div(&d, &g);
    KFrac ic = kfrac_inv(&dq.a[dq.n - 1], d.kind);       // monic denominator
    K2Poly icp = k2_const(ic, d.kind);
    K2Frac f = { k2_mul(&nq, &icp), k2_mul(&dq, &icp) };
    k2_free(&n); k2_free(&d); k2_free(&g); k2_free(&nq); k2_free(&dq); k2_free(&icp);
    return f;
}

// d = 1: the monic denominator of a polynomial
static bool k2_is_one(const K2Poly* d) { return d->n == 1 && kp_is_one(&d->a[0].n) && kp_is_one(&d->a[0].d); }

K2Frac k2frac_add(const K2Frac* a, const K2Frac* b, int sign) {
    if (k2_is_one(&a->d) && k2_is_one(&b->d)) return (K2Frac){ k2_add(&a->n, &b->n, sign), k2_copy(&a->d) }; // no gcd needed
    // Henrici: with g = gcd(a.d, b.d), t = a.n·(b.d/g) ± b.n·(a.d/g) can share factors only with g, so
    // the full gcd over K1 (whose coefficients swell) is replaced by gcd(t, g), which is 1 when g = 1
    K2Poly g = k2_gcd(&a->d, &b->d), ad = k2_exact_div(&a->d, &g), bd = k2_exact_div(&b->d, &g);
    K2Poly l = k2_mul(&a->n, &bd), r = k2_mul(&b->n, &ad), t = k2_add(&l, &r, sign);
    K2Frac f;
    if (t.n == 0) f = (K2Frac){ t, k2_const(kfrac_make(kp_one(), kp_one(), t.kind), t.kind) };
    else {
        K2Poly g2 = k2_gcd(&t, &g), bq = k2_exact_div(&b->d, &g2);
        f = (K2Frac){ k2_exact_div(&t, &g2), k2_mul(&ad, &bq) };   // monic: a product of monic quotients
        k2_free(&g2); k2_free(&bq); k2_free(&t);
    }
    k2_free(&g); k2_free(&ad); k2_free(&bd); k2_free(&l); k2_free(&r);
    return f;
}

K2Frac k2frac_mul(const K2Frac* a, const K2Frac* b) {
    if (k2_is_one(&a->d) && k2_is_one(&b->d)) return (K2Frac){ k2_mul(&a->n, &b->n), k2_copy(&a->d) };
    return k2frac_make(k2_mul(&a->n, &b->n), k2_mul(&a->d, &b->d));
}

K2Frac k2frac_inv(const K2Frac* a) { return k2frac_make(k2_copy(&a->d), k2_copy(&a->n)); }

void k2frac_free(K2Frac* f) { k2_free(&f->n); k2_free(&f->d); }
