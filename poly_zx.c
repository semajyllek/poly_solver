#include "poly_zx.h"
#include <stdlib.h>

// --- lifecycle ---

ZxPoly zx_alloc(int deg) {
    ZxPoly a;
    a.deg = deg;
    int n = deg < 0 ? 1 : deg + 1;
    a.c = malloc((size_t)n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) mpz_init(a.c[i]);
    return a;
}

// invariant: exactly max(deg,0)+1 coefficients are initialized; zx_trim clears what it drops

ZxPoly zx_copy(const ZxPoly* a) {
    ZxPoly r = zx_alloc(a->deg);
    for (int i = 0; i <= a->deg; i++) mpz_set(r.c[i], a->c[i]);
    return r;
}

void zx_free(ZxPoly* a) {
    int n = a->deg < 0 ? 1 : a->deg + 1;
    for (int i = 0; i < n; i++) mpz_clear(a->c[i]);
    free(a->c);
    a->c = NULL;
    a->deg = -1;
}

void zx_trim(ZxPoly* a) {
    // clear dropped leading zeros so zx_free's count (deg+1) stays exact
    while (a->deg > 0 && mpz_sgn(a->c[a->deg]) == 0) {
        mpz_clear(a->c[a->deg]);
        a->deg--;
    }
    if (a->deg == 0 && mpz_sgn(a->c[0]) == 0) a->deg = -1; // keep c[0] allocated
}

ZxPoly zx_from_poly(const Polynomial* p) {
    ZxPoly a = zx_alloc(polynomial_degree(p));
    for (Term* t = p->head; t; t = t->next)
        mpz_set(a.c[t->exponent], mpq_numref(t->coeff.val));
    return a;
}

Polynomial* zx_to_poly(const ZxPoly* a) {
    Polynomial* p = create_polynomial();
    Rational r;
    rat_init(&r);
    for (int i = a->deg; i >= 0; i--) {
        mpq_set_z(r.val, a->c[i]);
        add_term(p, &r, i);
    }
    rat_clear(&r);
    return p;
}

// --- arithmetic ---

ZxPoly zx_mul(const ZxPoly* a, const ZxPoly* b) {
    if (a->deg < 0 || b->deg < 0) return zx_alloc(-1);
    ZxPoly r = zx_alloc(a->deg + b->deg);
    for (int i = 0; i <= a->deg; i++) {
        if (mpz_sgn(a->c[i]) == 0) continue;
        for (int j = 0; j <= b->deg; j++) mpz_addmul(r.c[i + j], a->c[i], b->c[j]);
    }
    return r;
}

bool zx_divexact(const ZxPoly* a, const ZxPoly* b, ZxPoly* q) {
    if (b->deg < 0 || a->deg < b->deg) return false;
    ZxPoly r = zx_copy(a);
    *q = zx_alloc(a->deg - b->deg);
    mpz_t t;
    mpz_init(t);
    bool ok = true;
    for (int i = a->deg - b->deg; i >= 0 && ok; i--) {
        // leading coefficient must divide exactly, otherwise b does not divide a
        if (!mpz_divisible_p(r.c[b->deg + i], b->c[b->deg])) { ok = false; break; }
        mpz_divexact(q->c[i], r.c[b->deg + i], b->c[b->deg]);
        for (int j = 0; j <= b->deg; j++) {
            mpz_mul(t, q->c[i], b->c[j]);
            mpz_sub(r.c[i + j], r.c[i + j], t);
        }
    }
    for (int i = 0; i < b->deg && ok; i++) // remainder must vanish
        if (mpz_sgn(r.c[i]) != 0) ok = false;
    mpz_clear(t);
    zx_free(&r);
    if (!ok) zx_free(q);
    return ok;
}

void zx_make_primitive(ZxPoly* a) {
    if (a->deg < 0) return;
    mpz_t g;
    mpz_init(g);
    for (int i = 0; i <= a->deg; i++) mpz_gcd(g, g, a->c[i]);
    if (mpz_sgn(a->c[a->deg]) < 0) mpz_neg(g, g); // also normalizes the sign
    for (int i = 0; i <= a->deg; i++) mpz_divexact(a->c[i], a->c[i], g);
    mpz_clear(g);
}

void zx_mods(ZxPoly* a, const mpz_t m) {
    mpz_t half;
    mpz_init(half);
    mpz_fdiv_q_2exp(half, m, 1);
    for (int i = 0; i <= a->deg; i++) {
        mpz_mod(a->c[i], a->c[i], m);
        if (mpz_cmp(a->c[i], half) > 0) mpz_sub(a->c[i], a->c[i], m);
    }
    mpz_clear(half);
    zx_trim(a);
}

ZpPoly zx_to_zp(const ZxPoly* a, const ZpField* F) {
    ZpPoly r = zp_alloc(a->deg);
    for (int i = 0; i <= a->deg; i++) r.c[i] = mpz_fdiv_ui(a->c[i], F->p);
    zp_trim(&r);
    return r;
}

ZxPoly zx_from_zp(const ZpPoly* a) {
    ZxPoly r = zx_alloc(a->deg);
    for (int i = 0; i <= a->deg; i++) mpz_set_ui(r.c[i], a->c[i]);
    return r;
}

void zx_factor_list_free(ZxFactorList* L) {
    for (int i = 0; i < L->count; i++) zx_free(&L->items[i]);
    free(L->items);
    L->items = NULL;
    L->count = 0;
}
