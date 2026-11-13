#include "roots.h"
#include "poly_zx.h"
#include "factor.h"
#include <stdlib.h>
#include <complex.h>

// t-interval (c/2^k, (c+1)/2^k) inside (0,1), the unit interval of a scaled polynomial
typedef struct {
    mpz_t c;
    unsigned long k;
} Dyadic;

// maps t in (0,1) back to x = sign · B · t and collects results
typedef struct {
    const Rational* B;
    int sign;
    RootIntervalList* out;
} IsoCtx;

static void push_interval(RootIntervalList* L, const Rational* lo, const Rational* hi) {
    L->items = realloc(L->items, (size_t)(L->count + 1) * sizeof(RootInterval));
    rat_init_set(&L->items[L->count].lo, lo);
    rat_init_set(&L->items[L->count].hi, hi);
    L->count++;
}

// --- integer polynomial transforms ---

// g(x) -> g(x + 1), in place (repeated synthetic division by x - 1)
static void taylor_shift1(ZxPoly* g) {
    for (int i = 0; i < g->deg; i++)
        for (int j = g->deg - 1; j >= i; j--) mpz_add(g->c[j], g->c[j], g->c[j + 1]);
}

// g(x) -> 2^n · g(x/2), keeping integer coefficients
static ZxPoly scale_half(const ZxPoly* g) {
    ZxPoly r = zx_copy(g);
    for (int i = 0; i <= r.deg; i++) mpz_mul_2exp(r.c[i], r.c[i], (mp_bitcnt_t)(r.deg - i));
    return r;
}

// Descartes bound on the number of roots of g in (0,1):
// sign variations of (x+1)^n · g(1/(x+1))
static int descartes_unit(const ZxPoly* g) {
    ZxPoly h = zx_alloc(g->deg);
    for (int i = 0; i <= g->deg; i++) mpz_set(h.c[i], g->c[g->deg - i]); // reverse
    taylor_shift1(&h);
    int var = 0, last = 0;
    for (int i = 0; i <= h.deg; i++) {
        int s = mpz_sgn(h.c[i]);
        if (s != 0 && last != 0 && s != last) var++;
        if (s != 0) last = s;
    }
    zx_free(&h);
    return var;
}

// dyadic point c/2^k scaled into x-space: sign · B · c / 2^k
static void to_x(Rational* out, const Dyadic* d, const IsoCtx* ctx) {
    mpq_set_z(out->val, d->c);
    mpq_div_2exp(out->val, out->val, (mp_bitcnt_t)d->k);
    rat_mul(out, out, ctx->B);
    if (ctx->sign < 0) rat_neg(out, out);
}

// record the x-interval corresponding to the t-interval of d
static void record(const Dyadic* d, IsoCtx* ctx) {
    Rational a, b;
    rat_init(&a); rat_init(&b);
    Dyadic right;
    mpz_init(right.c);
    mpz_add_ui(right.c, d->c, 1);
    right.k = d->k;
    to_x(&a, d, ctx);
    to_x(&b, &right, ctx);
    if (rat_cmp(&a, &b) > 0) push_interval(ctx->out, &b, &a);
    else push_interval(ctx->out, &a, &b);
    mpz_clear(right.c);
    rat_clear(&a); rat_clear(&b);
}

// record an exact root at the dyadic point d
static void record_exact(const Dyadic* d, IsoCtx* ctx) {
    Rational x;
    rat_init(&x);
    to_x(&x, d, ctx);
    push_interval(ctx->out, &x, &x);
    rat_clear(&x);
}

// divide out a root at x = 0 (shift coefficients down one place)
static void drop_zero_root(ZxPoly* g) {
    for (int i = 0; i < g->deg; i++) mpz_set(g->c[i], g->c[i + 1]);
    mpz_set_ui(g->c[g->deg], 0);
    zx_trim(g);
}

// bisection driven by the Descartes count on the unit interval of g
static void isolate_unit(const ZxPoly* g, const Dyadic* d, IsoCtx* ctx) {
    int v = descartes_unit(g);
    if (v == 0) return;
    if (v == 1) { record(d, ctx); return; }
    ZxPoly left = scale_half(g);          // covers (0, 1/2)
    ZxPoly right = zx_copy(&left);
    taylor_shift1(&right);                // covers (1/2, 1)
    Dyadic dl, dr;
    mpz_init(dl.c); mpz_init(dr.c);
    mpz_mul_2exp(dl.c, d->c, 1);
    mpz_add_ui(dr.c, dl.c, 1);
    dl.k = dr.k = d->k + 1;
    if (right.deg >= 0 && mpz_sgn(right.c[0]) == 0) { // the midpoint is an exact root
        record_exact(&dr, ctx);
        drop_zero_root(&right);
    }
    isolate_unit(&left, &dl, ctx);
    isolate_unit(&right, &dr, ctx);
    mpz_clear(dl.c); mpz_clear(dr.c);
    zx_free(&left); zx_free(&right);
}

// power of two strictly above every root magnitude (Cauchy bound 1 + max|a_i/a_n|)
static void root_bound(Rational* B, const ZxPoly* f) {
    mpz_t mx, q;
    mpz_init(mx); mpz_init(q);
    for (int i = 0; i < f->deg; i++)
        if (mpz_cmpabs(f->c[i], mx) > 0) mpz_abs(mx, f->c[i]);
    mpz_abs(q, f->c[f->deg]);
    mpz_cdiv_q(mx, mx, q);
    mpz_add_ui(mx, mx, 1);
    mpz_set_ui(q, 1);
    while (mpz_cmp(q, mx) <= 0) mpz_mul_2exp(q, q, 1);
    mpq_set_z(B->val, q);
    mpz_clear(mx); mpz_clear(q);
}

// g(t) = f(sign · B · t), so roots of f on one side of 0 land in (0,1)
static ZxPoly scale_to_unit(const ZxPoly* f, const Rational* B, int sign) {
    ZxPoly g = zx_copy(f);
    mpz_t pw;
    mpz_init_set_ui(pw, 1);
    for (int i = 0; i <= g.deg; i++) {
        mpz_mul(g.c[i], g.c[i], pw);
        mpz_mul(pw, pw, mpq_numref(B->val));
        if (sign < 0) mpz_neg(pw, pw);
    }
    mpz_clear(pw);
    return g;
}

static void sort_intervals(RootIntervalList* L) {
    for (int i = 1; i < L->count; i++)
        for (int j = i; j > 0 && rat_cmp(&L->items[j - 1].lo, &L->items[j].lo) > 0; j--) {
            RootInterval t = L->items[j];
            L->items[j] = L->items[j - 1];
            L->items[j - 1] = t;
        }
}

RootIntervalList isolate_real_roots(const Polynomial* f) {
    RootIntervalList L = { NULL, 0 };
    if (polynomial_degree(f) <= 0) return L;
    Polynomial* prim = poly_primitive_part(f);
    ZxPoly z = zx_from_poly(prim);
    free_polynomial(prim);
    Rational B;
    rat_init(&B);
    if (mpz_sgn(z.c[0]) == 0) { // root at exactly 0
        Rational zero = rat_zero_val();
        push_interval(&L, &zero, &zero);
        rat_clear(&zero);
        drop_zero_root(&z);
    }
    if (z.deg > 0) {
        root_bound(&B, &z);
        Dyadic unit;
        mpz_init(unit.c);
        unit.k = 0;
        for (int sign = 1; sign >= -1; sign -= 2) { // positive side, then negative side
            ZxPoly g = scale_to_unit(&z, &B, sign);
            IsoCtx ctx = { &B, sign, &L };
            isolate_unit(&g, &unit, &ctx);
            zx_free(&g);
        }
        mpz_clear(unit.c);
    }
    sort_intervals(&L);
    rat_clear(&B);
    zx_free(&z);
    return L;
}

void refine_root(RootInterval* r, const Polynomial* f, const Rational* width) {
    Rational v, mid, w;
    rat_init(&v); rat_init(&mid); rat_init(&w);
    poly_eval_rat(&v, f, &r->lo);
    int sign_lo = mpq_sgn(v.val);
    rat_sub(&w, &r->hi, &r->lo);
    // endpoints of an isolating interval are never roots, so f changes sign across it
    while (sign_lo != 0 && rat_cmp(&w, width) > 0) {
        rat_add(&mid, &r->lo, &r->hi);
        mpq_div_2exp(mid.val, mid.val, 1);
        poly_eval_rat(&v, f, &mid);
        int s = mpq_sgn(v.val);
        if (s == 0) { rat_set(&r->lo, &mid); rat_set(&r->hi, &mid); break; }
        if (s == sign_lo) rat_set(&r->lo, &mid);
        else rat_set(&r->hi, &mid);
        rat_sub(&w, &r->hi, &r->lo);
    }
    rat_clear(&v); rat_clear(&mid); rat_clear(&w);
}

void root_interval_list_free(RootIntervalList* L) {
    for (int i = 0; i < L->count; i++) {
        rat_clear(&L->items[i].lo);
        rat_clear(&L->items[i].hi);
    }
    free(L->items);
    L->items = NULL;
    L->count = 0;
}

// --- numerical complex roots (Durand–Kerner / Weierstrass iteration) ---

#define DK_ITERATIONS 2000

void approx_complex_roots(const Polynomial* p, double complex* out) {
    int n = polynomial_degree(p);
    double* c = malloc((size_t)(n + 1) * sizeof(double)); // monic coefficients, c[k] of x^k
    Rational r;
    rat_init(&r);
    for (int k = 0; k <= n; k++) { poly_coeff(&r, p, k); c[k] = rat_to_double(&r); }
    rat_clear(&r);
    for (int k = 0; k <= n; k++) c[k] /= c[n] == 0 ? 1 : c[n];
    for (int i = 0; i < n; i++) out[i] = cpow(0.4 + 0.9 * I, i); // standard distinct starts
    for (int it = 0; it < DK_ITERATIONS; it++) {
        double change = 0;
        for (int i = 0; i < n; i++) {
            double complex v = 1, den = 1;
            for (int k = n - 1; k >= 0; k--) v = v * out[i] + c[k]; // Horner, monic
            for (int j = 0; j < n; j++) if (j != i) den *= out[i] - out[j];
            double complex step = v / den;
            out[i] -= step;
            change += cabs(step);
        }
        if (change < 1e-15) break;
    }
    free(c);
}
