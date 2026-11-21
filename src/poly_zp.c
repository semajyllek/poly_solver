#include "poly_zp.h"
#include <stdlib.h>
#include <string.h>

// --- lifecycle ---

ZpPoly zp_alloc(int deg) {
    ZpPoly a;
    a.deg = deg;
    a.c = calloc((size_t)(deg < 0 ? 1 : deg + 1), sizeof(unsigned long));
    return a;
}

ZpPoly zp_copy(const ZpPoly* a) {
    ZpPoly r = zp_alloc(a->deg);
    if (a->deg >= 0) memcpy(r.c, a->c, (size_t)(a->deg + 1) * sizeof(unsigned long));
    return r;
}

void zp_free(ZpPoly* a) {
    free(a->c);
    a->c = NULL;
    a->deg = -1;
}

// drop leading zero coefficients
void zp_trim(ZpPoly* a) {
    while (a->deg >= 0 && a->c[a->deg] == 0) a->deg--;
}

// --- scalar arithmetic (p < 2^31, so products fit in 64 bits) ---

static unsigned long mulm(unsigned long a, unsigned long b, const ZpField* F) {
    return (unsigned long)(((unsigned long long)a * b) % F->p);
}

unsigned long zp_inv(unsigned long a, const ZpField* F) {
    // extended Euclid on integers
    long long t = 0, nt = 1, r = (long long)F->p, nr = (long long)(a % F->p);
    while (nr != 0) {
        long long q = r / nr, tmp;
        tmp = t - q * nt; t = nt; nt = tmp;
        tmp = r - q * nr; r = nr; nr = tmp;
    }
    if (t < 0) t += (long long)F->p;
    return (unsigned long)t;
}

// --- polynomial arithmetic ---

ZpPoly zp_add(const ZpPoly* a, const ZpPoly* b, const ZpField* F) {
    int d = a->deg > b->deg ? a->deg : b->deg;
    ZpPoly r = zp_alloc(d);
    for (int i = 0; i <= d; i++) {
        unsigned long x = i <= a->deg ? a->c[i] : 0, y = i <= b->deg ? b->c[i] : 0;
        r.c[i] = (x + y) % F->p;
    }
    zp_trim(&r);
    return r;
}

ZpPoly zp_sub(const ZpPoly* a, const ZpPoly* b, const ZpField* F) {
    int d = a->deg > b->deg ? a->deg : b->deg;
    ZpPoly r = zp_alloc(d);
    for (int i = 0; i <= d; i++) {
        unsigned long x = i <= a->deg ? a->c[i] : 0, y = i <= b->deg ? b->c[i] : 0;
        r.c[i] = (x + F->p - y) % F->p;
    }
    zp_trim(&r);
    return r;
}

ZpPoly zp_mul(const ZpPoly* a, const ZpPoly* b, const ZpField* F) {
    if (a->deg < 0 || b->deg < 0) return zp_alloc(-1);
    ZpPoly r = zp_alloc(a->deg + b->deg);
    for (int i = 0; i <= a->deg; i++) {
        if (a->c[i] == 0) continue;
        for (int j = 0; j <= b->deg; j++)
            r.c[i + j] = (r.c[i + j] + mulm(a->c[i], b->c[j], F)) % F->p;
    }
    zp_trim(&r);
    return r;
}

typedef struct {
    ZpPoly q, r;
} ZpDivResult;

// long division a = q·b + r
static ZpDivResult zp_divrem(const ZpPoly* a, const ZpPoly* b, const ZpField* F) {
    ZpPoly r = zp_copy(a);
    unsigned long inv_lc = zp_inv(b->c[b->deg], F);
    int qdeg = a->deg - b->deg;
    ZpPoly q = zp_alloc(qdeg < 0 ? -1 : qdeg);
    for (int i = qdeg; i >= 0; i--) {
        unsigned long coef = mulm(r.c[b->deg + i], inv_lc, F);
        q.c[i] = coef;
        if (coef == 0) continue;
        for (int j = 0; j <= b->deg; j++) // r -= coef * x^i * b
            r.c[i + j] = (r.c[i + j] + F->p - mulm(coef, b->c[j], F)) % F->p;
    }
    zp_trim(&r);
    zp_trim(&q);
    return (ZpDivResult){ q, r };
}

ZpPoly zp_rem(const ZpPoly* a, const ZpPoly* b, const ZpField* F) {
    ZpDivResult d = zp_divrem(a, b, F);
    zp_free(&d.q);
    return d.r;
}

ZpPoly zp_quo(const ZpPoly* a, const ZpPoly* b, const ZpField* F) {
    ZpDivResult d = zp_divrem(a, b, F);
    zp_free(&d.r);
    return d.q;
}

void zp_make_monic(ZpPoly* a, const ZpField* F) {
    if (a->deg < 0) return;
    unsigned long inv = zp_inv(a->c[a->deg], F);
    for (int i = 0; i <= a->deg; i++) a->c[i] = mulm(a->c[i], inv, F);
}

ZpPoly zp_gcd(const ZpPoly* a, const ZpPoly* b, const ZpField* F) {
    ZpPoly x = zp_copy(a), y = zp_copy(b);
    while (y.deg >= 0) {
        ZpPoly r = zp_rem(&x, &y, F);
        zp_free(&x);
        x = y;
        y = r;
    }
    zp_free(&y);
    zp_make_monic(&x, F);
    return x;
}

ZpPoly zp_deriv(const ZpPoly* a, const ZpField* F) {
    if (a->deg <= 0) return zp_alloc(-1);
    ZpPoly r = zp_alloc(a->deg - 1);
    for (int i = 1; i <= a->deg; i++) r.c[i - 1] = mulm(a->c[i], (unsigned long)i % F->p, F);
    zp_trim(&r);
    return r;
}

// (a * b) mod R->mod
static ZpPoly zp_mulmod(const ZpPoly* a, const ZpPoly* b, const ZpRing* R) {
    ZpPoly prod = zp_mul(a, b, R->F);
    ZpPoly r = zp_rem(&prod, R->mod, R->F);
    zp_free(&prod);
    return r;
}

ZpPoly zp_powmod(const ZpPoly* base, const mpz_t e, const ZpRing* R) {
    ZpPoly result = zp_alloc(0);
    result.c[0] = 1;
    ZpPoly b = zp_rem(base, R->mod, R->F);
    // left-to-right square and multiply over the bits of e
    for (long i = (long)mpz_sizeinbase(e, 2) - 1; i >= 0 && mpz_sgn(e) > 0; i--) {
        ZpPoly sq = zp_mulmod(&result, &result, R);
        zp_free(&result);
        result = sq;
        if (mpz_tstbit(e, (mp_bitcnt_t)i)) {
            ZpPoly m = zp_mulmod(&result, &b, R);
            zp_free(&result);
            result = m;
        }
    }
    zp_free(&b);
    return result;
}

// scale polynomial by a constant
static void zp_scale(ZpPoly* a, unsigned long k, const ZpField* F) {
    for (int i = 0; i <= a->deg; i++) a->c[i] = mulm(a->c[i], k, F);
    zp_trim(a);
}

typedef struct {
    const ZpPoly* q; // quotient of the current division step
    const ZpField* F;
} ZpRowCtx;

// old - q·cur (one row of the extended Euclidean recurrence)
static ZpPoly zp_euclid_row(const ZpPoly* old, const ZpPoly* cur, const ZpRowCtx* ctx) {
    ZpPoly qc = zp_mul(ctx->q, cur, ctx->F);
    ZpPoly r = zp_sub(old, &qc, ctx->F);
    zp_free(&qc);
    return r;
}

ZpBezout zp_xgcd(const ZpPoly* a, const ZpPoly* b, const ZpField* F) {
    ZpPoly r0 = zp_copy(a), r1 = zp_copy(b);
    ZpPoly s0 = zp_alloc(0), s1 = zp_alloc(-1);
    ZpPoly t0 = zp_alloc(-1), t1 = zp_alloc(0);
    s0.c[0] = 1;
    t1.c[0] = 1;
    while (r1.deg >= 0) {
        ZpDivResult d = zp_divrem(&r0, &r1, F);
        ZpRowCtx row = { &d.q, F };
        ZpPoly s2 = zp_euclid_row(&s0, &s1, &row);
        ZpPoly t2 = zp_euclid_row(&t0, &t1, &row);
        zp_free(&d.q); zp_free(&r0); zp_free(&s0); zp_free(&t0);
        r0 = r1; r1 = d.r;
        s0 = s1; s1 = s2;
        t0 = t1; t1 = t2;
    }
    // r0 is a nonzero constant when gcd(a,b)=1: normalize so s*a + t*b = 1
    unsigned long inv = zp_inv(r0.c[0], F);
    zp_scale(&s0, inv, F);
    zp_scale(&t0, inv, F);
    zp_free(&r0); zp_free(&r1); zp_free(&s1); zp_free(&t1);
    return (ZpBezout){ s0, t0 };
}

// --- factoring over F_p ---

typedef struct {
    const ZpField* F;
    unsigned long long rng; // deterministic LCG state, so output is reproducible
    ZpFactorList* out;
} ZpFactorCtx;

static void list_push(ZpFactorList* L, ZpPoly f) {
    L->items = realloc(L->items, (size_t)(L->count + 1) * sizeof(ZpPoly));
    L->items[L->count++] = f;
}

// random polynomial of degree < n
static ZpPoly random_poly(int n, ZpFactorCtx* ctx) {
    ZpPoly a = zp_alloc(n - 1);
    for (int i = 0; i < n; i++) {
        ctx->rng = ctx->rng * 6364136223846793005ULL + 1442695040888963407ULL;
        a.c[i] = (unsigned long)((ctx->rng >> 33) % ctx->F->p);
    }
    zp_trim(&a);
    return a;
}

// one Cantor–Zassenhaus attempt: gcd(a^((p^d-1)/2) - 1, g)
static ZpPoly cz_try_split(const ZpPoly* g, int d, ZpFactorCtx* ctx) {
    mpz_t e;
    mpz_init(e);
    mpz_ui_pow_ui(e, ctx->F->p, (unsigned long)d);
    mpz_sub_ui(e, e, 1);
    mpz_fdiv_q_2exp(e, e, 1);
    ZpPoly a = random_poly(g->deg, ctx);
    ZpRing R = { ctx->F, g };
    ZpPoly b = zp_powmod(&a, e, &R);
    if (b.deg < 0) { zp_free(&b); b = zp_alloc(0); }
    b.c[0] = (b.c[0] + ctx->F->p - 1) % ctx->F->p; // b - 1
    zp_trim(&b);
    ZpPoly u = zp_gcd(&b, g, ctx->F);
    zp_free(&a); zp_free(&b);
    mpz_clear(e);
    return u;
}

// split g (product of irreducibles all of degree d) into its factors
static void equal_degree(const ZpPoly* g, int d, ZpFactorCtx* ctx) {
    if (g->deg == d) { list_push(ctx->out, zp_copy(g)); return; }
    for (;;) {
        ZpPoly u = cz_try_split(g, d, ctx);
        if (u.deg > 0 && u.deg < g->deg) {
            ZpPoly v = zp_quo(g, &u, ctx->F);
            equal_degree(&u, d, ctx);
            equal_degree(&v, d, ctx);
            zp_free(&u); zp_free(&v);
            return;
        }
        zp_free(&u);
    }
}

// x as a polynomial
static ZpPoly zp_x(void) {
    ZpPoly x = zp_alloc(1);
    x.c[1] = 1;
    return x;
}

ZpFactorList zp_factor_squarefree(const ZpPoly* f, const ZpField* F) {
    ZpFactorList L = { NULL, 0 };
    ZpFactorCtx ctx = { F, 0x2545F4914F6CDD1DULL, &L };
    ZpPoly rest = zp_copy(f), x = zp_x(), h = zp_x();
    zp_make_monic(&rest, F);
    mpz_t p;
    mpz_init_set_ui(p, F->p);
    // distinct-degree: gcd(x^(p^d) - x, rest) collects the degree-d factors
    for (int d = 1; 2 * d <= rest.deg; d++) {
        ZpRing R = { F, &rest };
        ZpPoly hp = zp_powmod(&h, p, &R);
        zp_free(&h);
        h = hp;
        ZpPoly hx = zp_sub(&h, &x, F);
        ZpPoly g = zp_gcd(&hx, &rest, F);
        zp_free(&hx);
        if (g.deg > 0) {
            equal_degree(&g, d, &ctx);
            ZpPoly q = zp_quo(&rest, &g, F);
            zp_free(&rest);
            rest = q;
            ZpPoly hr = zp_rem(&h, &rest, F);
            zp_free(&h);
            h = hr;
        }
        zp_free(&g);
    }
    if (rest.deg > 0) list_push(&L, zp_copy(&rest)); // what remains is irreducible
    zp_free(&rest); zp_free(&x); zp_free(&h);
    mpz_clear(p);
    return L;
}

void zp_factor_list_free(ZpFactorList* L) {
    for (int i = 0; i < L->count; i++) zp_free(&L->items[i]);
    free(L->items);
    L->items = NULL;
    L->count = 0;
}
