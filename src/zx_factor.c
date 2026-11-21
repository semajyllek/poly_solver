// Berlekamp–Zassenhaus factorization over Z[x]:
//   1. pick a small prime p with f square-free mod p (fewest modular factors of a few tries)
//   2. factor f mod p (Cantor–Zassenhaus, poly_zp.c)
//   3. Hensel-lift the monic factors to p^k > 2·|lc|·(Mignotte bound)
//   4. recombine subsets; a candidate is accepted only if it divides f exactly in Z[x],
//      so a lifting bug can only cause under-factoring, never a wrong factor
#include "poly_zx.h"
#include <stdlib.h>

#define PRIME_TRIES 5 // valid primes examined before settling on the best one

// --- prime selection ---

static bool is_prime_ul(unsigned long n) {
    if (n < 2) return false;
    for (unsigned long d = 2; d * d <= n; d++)
        if (n % d == 0) return false;
    return true;
}

// true when p does not divide lc(f) and f stays square-free mod p
static bool prime_is_usable(const ZxPoly* f, const ZpField* F) {
    if (mpz_fdiv_ui(f->c[f->deg], F->p) == 0) return false;
    ZpPoly fp = zx_to_zp(f, F);
    ZpPoly d = zp_deriv(&fp, F);
    ZpPoly g = zp_gcd(&fp, &d, F);
    bool ok = (g.deg == 0);
    zp_free(&fp); zp_free(&d); zp_free(&g);
    return ok;
}

// factor f mod several primes and keep the one with the fewest factors
static ZpFactorList choose_prime(const ZxPoly* f, ZpField* best) {
    ZpFactorList best_list = { NULL, 0 };
    int tried = 0;
    for (unsigned long p = 3; tried < PRIME_TRIES; p += 2) {
        if (!is_prime_ul(p)) continue;
        ZpField F = { p };
        if (!prime_is_usable(f, &F)) continue;
        tried++;
        ZpPoly fp = zx_to_zp(f, &F);
        ZpFactorList L = zp_factor_squarefree(&fp, &F);
        zp_free(&fp);
        if (best_list.items == NULL || L.count < best_list.count) {
            zp_factor_list_free(&best_list);
            best_list = L;
            *best = F;
        } else {
            zp_factor_list_free(&L);
        }
        if (best_list.count == 1) break; // irreducible mod p => irreducible over Z
    }
    return best_list;
}

// --- coefficient bound ---

// smallest k with p^k > 2·|lc|·(n+1)·2^n·max|a_i|  (Mignotte, loosened by sqrt(n+1) <= n+1)
static void lifting_modulus(mpz_t P, const ZxPoly* f, unsigned long p) {
    mpz_t bound, mx;
    mpz_init(bound);
    mpz_init(mx);
    for (int i = 0; i <= f->deg; i++)
        if (mpz_cmpabs(f->c[i], mx) > 0) mpz_abs(mx, f->c[i]);
    mpz_mul_2exp(bound, mx, (mp_bitcnt_t)f->deg + 1); // 2·2^n·max
    mpz_mul_ui(bound, bound, (unsigned long)f->deg + 1);
    mpz_mul(bound, bound, f->c[f->deg]);
    mpz_abs(bound, bound);
    mpz_set_ui(P, p);
    while (mpz_cmp(P, bound) <= 0) mpz_mul_ui(P, P, p);
    mpz_clear(bound);
    mpz_clear(mx);
}

// --- Hensel lifting ---

typedef struct {
    ZxPoly g, h;   // monic lifts with f ≡ g·h (mod M)
    ZpPoly s, t;   // s·g + t·h ≡ 1 (mod p)
    mpz_t M;       // current modulus p^j
} HenselPair;

typedef struct {
    const ZpField* F;
    mpz_srcptr P;  // final modulus p^k
    ZxPoly* out;   // lifted factors, filled in order
    int next;      // next free slot in out
} LiftCtx;

// e = ((target - g·h) mod M·p) / M, reduced mod p
static ZpPoly hensel_error(const HenselPair* hp, const ZxPoly* target, const ZpField* F) {
    ZxPoly gh = zx_mul(&hp->g, &hp->h);
    int d = target->deg > gh.deg ? target->deg : gh.deg;
    ZxPoly diff = zx_alloc(d);
    for (int i = 0; i <= d; i++) {
        if (i <= target->deg) mpz_add(diff.c[i], diff.c[i], target->c[i]);
        if (i <= gh.deg) mpz_sub(diff.c[i], diff.c[i], gh.c[i]);
        mpz_fdiv_q(diff.c[i], diff.c[i], hp->M); // exact up to the residue we keep mod p
    }
    zx_trim(&diff);
    ZpPoly e = zx_to_zp(&diff, F);
    zx_free(&gh);
    zx_free(&diff);
    return e;
}

// add M·delta to a lifted factor
static void add_scaled(ZxPoly* a, const ZpPoly* delta, const mpz_t M) {
    for (int i = 0; i <= delta->deg && i <= a->deg; i++)
        mpz_addmul_ui(a->c[i], M, delta->c[i]);
}

// one linear lifting step: modulus M -> M·p, keeping g and h monic
static void hensel_step(HenselPair* hp, const ZxPoly* target, const ZpField* F) {
    ZpPoly e = hensel_error(hp, target, F);
    ZpPoly gp = zx_to_zp(&hp->g, F), hpp = zx_to_zp(&hp->h, F);
    ZpPoly es = zp_mul(&e, &hp->s, F), et = zp_mul(&e, &hp->t, F);
    ZpPoly dh = zp_rem(&es, &hpp, F); // deg < deg h keeps h monic
    ZpPoly dg = zp_rem(&et, &gp, F);  // deg < deg g keeps g monic
    add_scaled(&hp->g, &dg, hp->M);
    add_scaled(&hp->h, &dh, hp->M);
    mpz_mul_ui(hp->M, hp->M, F->p);
    zp_free(&e); zp_free(&gp); zp_free(&hpp);
    zp_free(&es); zp_free(&et); zp_free(&dh); zp_free(&dg);
}

// product of a run of modular factors
static ZpPoly zp_product(const ZpPoly* items, int n, const ZpField* F) {
    ZpPoly acc = zp_alloc(0);
    acc.c[0] = 1;
    for (int i = 0; i < n; i++) {
        ZpPoly t = zp_mul(&acc, &items[i], F);
        zp_free(&acc);
        acc = t;
    }
    return acc;
}

typedef struct {
    const ZpPoly* items;
    int n;
} ZpSlice;

// lift target ≡ Π slice (mod p) to mod P by recursive binary splitting
static void lift_tree(const ZxPoly* target, ZpSlice slice, LiftCtx* ctx) {
    const ZpPoly* items = slice.items;
    int n = slice.n;
    if (n == 1) {
        ctx->out[ctx->next] = zx_copy(target);
        zx_mods(&ctx->out[ctx->next], ctx->P);
        ctx->next++;
        return;
    }
    int half = n / 2;
    ZpPoly gp = zp_product(items, half, ctx->F);
    ZpPoly hpp = zp_product(items + half, n - half, ctx->F);
    HenselPair hp;
    ZpBezout bz = zp_xgcd(&gp, &hpp, ctx->F);
    hp.g = zx_from_zp(&gp);
    hp.h = zx_from_zp(&hpp);
    hp.s = bz.s;
    hp.t = bz.t;
    mpz_init_set_ui(hp.M, ctx->F->p);
    while (mpz_cmp(hp.M, ctx->P) < 0) hensel_step(&hp, target, ctx->F);
    zx_mods(&hp.g, ctx->P);
    zx_mods(&hp.h, ctx->P);
    lift_tree(&hp.g, (ZpSlice){ items, half }, ctx);
    lift_tree(&hp.h, (ZpSlice){ items + half, n - half }, ctx);
    zp_free(&gp); zp_free(&hpp); zp_free(&hp.s); zp_free(&hp.t);
    zx_free(&hp.g); zx_free(&hp.h);
    mpz_clear(hp.M);
}

// f / lc(f) mod P: the monic polynomial that the modular factors multiply to
static ZxPoly monic_target(const ZxPoly* f, const mpz_t P) {
    ZxPoly t = zx_copy(f);
    mpz_t inv;
    mpz_init(inv);
    mpz_invert(inv, f->c[f->deg], P);
    for (int i = 0; i <= t.deg; i++) mpz_mul(t.c[i], t.c[i], inv);
    zx_mods(&t, P);
    mpz_set_ui(t.c[t.deg], 1); // exact 1, never -(P-1)
    mpz_clear(inv);
    return t;
}

// --- recombination ---

typedef struct {
    ZxPoly f;          // part of the input not yet split off
    ZxPoly* lifted;    // lifted modular factors
    bool* used;        // consumed by an accepted true factor
    int r;             // number of lifted factors
    mpz_srcptr P;
    ZxFactorList* out;
} Recombine;

static void push_factor(ZxFactorList* L, ZxPoly g) {
    L->items = realloc(L->items, (size_t)(L->count + 1) * sizeof(ZxPoly));
    L->items[L->count++] = g;
}

// lc(f) · Π lifted[idx] mod P, then primitive part
static ZxPoly candidate(const Recombine* R, const int* idx, int s) {
    ZxPoly g = zx_alloc(0);
    mpz_set(g.c[0], R->f.c[R->f.deg]);
    for (int i = 0; i < s; i++) {
        ZxPoly t = zx_mul(&g, &R->lifted[idx[i]]);
        zx_free(&g);
        g = t;
        zx_mods(&g, R->P);
    }
    zx_make_primitive(&g);
    return g;
}

// next s-combination of {0..m-1} in lexicographic order; false when exhausted
static bool next_combination(int* c, int s, int m) {
    int i = s - 1;
    while (i >= 0 && c[i] == m - s + i) i--;
    if (i < 0) return false;
    c[i]++;
    for (int j = i + 1; j < s; j++) c[j] = c[j - 1] + 1;
    return true;
}

// indices of the lifted factors not yet consumed; returns how many
static int unused_indices(const Recombine* R, int* avail) {
    int m = 0;
    for (int i = 0; i < R->r; i++) if (!R->used[i]) avail[m++] = i;
    return m;
}

// try every s-subset of unused factors; on a true factor, split it off and return true
static bool try_size(Recombine* R, int s) {
    int* avail = malloc((size_t)R->r * sizeof(int));
    int* c = malloc((size_t)s * sizeof(int));
    int* idx = malloc((size_t)s * sizeof(int));
    int m = unused_indices(R, avail);
    for (int i = 0; i < s; i++) c[i] = i;
    bool found = false;
    do {
        for (int i = 0; i < s; i++) idx[i] = avail[c[i]];
        ZxPoly g = candidate(R, idx, s);
        ZxPoly q;
        if (g.deg > 0 && zx_divexact(&R->f, &g, &q)) {
            push_factor(R->out, g);
            zx_free(&R->f);
            R->f = q;
            for (int i = 0; i < s; i++) R->used[idx[i]] = true;
            found = true;
        } else {
            zx_free(&g);
        }
    } while (!found && next_combination(c, s, m));
    free(avail); free(c); free(idx);
    return found;
}

static int unused_count(const Recombine* R) {
    int n = 0;
    for (int i = 0; i < R->r; i++) n += !R->used[i];
    return n;
}

static void recombine(Recombine* R) {
    for (int s = 1; 2 * s <= unused_count(R); ) {
        if (!try_size(R, s)) s++; // same size again after a success
    }
    if (R->f.deg > 0) push_factor(R->out, zx_copy(&R->f)); // the rest is irreducible
}

// --- driver ---

ZxFactorList zx_factor_squarefree(const ZxPoly* f) {
    ZxFactorList out = { NULL, 0 };
    if (f->deg <= 1) { push_factor(&out, zx_copy(f)); return out; }
    ZpField F;
    ZpFactorList modular = choose_prime(f, &F);
    if (modular.count <= 1) {
        push_factor(&out, zx_copy(f));
        zp_factor_list_free(&modular);
        return out;
    }
    mpz_t P;
    mpz_init(P);
    lifting_modulus(P, f, F.p);
    ZxPoly target = monic_target(f, P);
    ZxPoly* lifted = malloc((size_t)modular.count * sizeof(ZxPoly));
    LiftCtx ctx = { &F, P, lifted, 0 };
    lift_tree(&target, (ZpSlice){ modular.items, modular.count }, &ctx);
    Recombine R = { zx_copy(f), lifted, calloc((size_t)modular.count, sizeof(bool)),
                    modular.count, P, &out };
    recombine(&R);
    for (int i = 0; i < modular.count; i++) zx_free(&lifted[i]);
    free(lifted);
    free(R.used);
    zx_free(&R.f);
    zx_free(&target);
    zp_factor_list_free(&modular);
    mpz_clear(P);
    return out;
}
