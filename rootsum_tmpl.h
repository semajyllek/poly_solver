// The log part for residues that are the roots α of an irreducible f ∈ Q[t] of degree n >= 3:
// Σ_{f(α)=0} α·ln(g_α), g_α = gcd(d, h - α·D(d)) over K(α)[θ] (Trager / Bronstein §5.6), written
// as a RootSum. Instantiated like logpair_tmpl.h (residue.c over K = Q(x), risch2.c over K1),
// with the same LT_ macros plus LT_C_COPY(x, a) and LT_C_IS_ZERO(a); include it after
// logpair_tmpl.h (it uses HD) and undefine the macros with logterm_undef.h.
//   K(α) = K[t]/(f): an element is an LT_P polynomial in t of degree < n (f stays irreducible
//   over K since its coefficients are constants), with the inverse from Bezout against f.
//   K(α)[θ]: kpoly_tmpl.h instantiated over those elements, for the monic gcd.
// g_α is found against e = gcd(d, Σ_j f_j·h^j·D(d)^(n-j) mod d), the product of the conjugate
// gcds, computed over K; e/n is the degree of each g_α. The derivative of the whole RootSum is
//   Σ_α α·D(g_α)/g_α = Tr(α·D(g_α)·(e/g_α)) / e,
// the trace Σ_i c_i·Tr(α^i) taken coefficientwise from the power sums of f's roots: an exact
// element of K(θ) for the caller's check. Generated: root_log, root_log_free.

// --- K(α) ---

typedef struct {
    LT_CTX x;
    LT_P F;      // f made monic, a polynomial in t over K
    int n;
    const Rational* ps; // ps[i] = Σ α^i over the conjugates, i < n (for traces)
} LT_PFX(Alg);

// p mod F (takes p)
static LT_P LT_PFX(a_reduce)(const LT_PFX(Alg)* A, LT_P* p) {
    LT_CTX x = A->x;
    (void)x;
    if (LT_P_DEG(p) < A->n) return *p;
    LT_KDIV r = LT_P_DIVMOD(x, p, &A->F);
    LT_P_FREE(p);
    LT_P_FREE(&r.q);
    return r.r;
}

static LT_P LT_PFX(a_mul)(const LT_PFX(Alg)* A, const LT_P* a, const LT_P* b) {
    LT_P p = LT_P_MUL(A->x, a, b);
    return LT_PFX(a_reduce)(A, &p);
}

// 1/a for a ≠ 0: s with s·a + t·F = 1
static LT_P LT_PFX(a_inv)(const LT_PFX(Alg)* A, const LT_P* a) {
    LT_CTX x = A->x;
    (void)x;
    LT_P one = LT_P_ONE(x);
    LT_KBEZ z = LT_P_BEZOUT(x, a, &A->F, &one);
    LT_P s = LT_P_COPY(x, &z.s);
    LT_KBEZ_FREE(&z);
    LT_P_FREE(&one);
    return s;
}

// k·a
static LT_P LT_PFX(a_scale_int)(const LT_PFX(Alg)* A, const LT_P* a, int k) {
    LT_CTX x = A->x;
    (void)x;
    Rational r = rat_from_int(k);
    LT_C c = LT_C_RAT(x, &r);
    LT_P s = LT_P_SCALE(x, a, &c);
    LT_C_FREE(&c);
    rat_clear(&r);
    return s;
}

// --- K(α)[θ]: Σ a[k]·θ^k, trimmed ---

typedef struct {
    LT_P* a;
    int n;
} LT_PFX(AP);

typedef struct {
    LT_PFX(AP) q, r;
} LT_PFX(ADiv);

typedef struct {
    LT_PFX(AP) s, t;
} LT_PFX(ABez);

typedef struct {
    LT_PFX(AP)* f;
    int m;
} LT_PFX(ASf);

static void LT_PFX(ap_free)(LT_PFX(AP)* p) {
    for (int i = 0; i < p->n; i++) LT_P_FREE(&p->a[i]);
    free(p->a);
    p->a = NULL; p->n = 0;
}

// Σ arr[i]·θ^i, taking the array; zero top coefficients dropped
static LT_PFX(AP) LT_PFX(ap_take)(LT_P* arr, int n) {
    LT_PFX(AP) p = { arr, n };
    while (p.n > 0 && LT_P_DEG(&p.a[p.n - 1]) < 0) LT_P_FREE(&p.a[--p.n]);
    return p;
}

static LT_PFX(AP) LT_PFX(ap_zeros)(const LT_PFX(Alg)* A, int n) {
    (void)A;
    LT_P* a = malloc((size_t)(n > 0 ? n : 1) * sizeof(LT_P));
    for (int i = 0; i < n; i++) a[i] = LT_P_ZERO(A->x);
    return (LT_PFX(AP)){ a, n };
}

// c·θ^k (c copied)
static LT_PFX(AP) LT_PFX(ap_term)(const LT_PFX(Alg)* A, const LT_P* c, int k) {
    LT_PFX(AP) p = LT_PFX(ap_zeros)(A, k + 1);
    LT_P_FREE(&p.a[k]);
    p.a[k] = LT_P_COPY(A->x, c);
    return LT_PFX(ap_take)(p.a, p.n);
}

static LT_PFX(AP) LT_PFX(ap_add)(const LT_PFX(Alg)* A, const LT_PFX(AP)* a, const LT_PFX(AP)* b, int sign) {
    int n = a->n > b->n ? a->n : b->n;
    LT_PFX(AP) r = LT_PFX(ap_zeros)(A, n);
    for (int i = 0; i < n; i++) {
        LT_P z = LT_P_ZERO(A->x);
        const LT_P *u = i < a->n ? &a->a[i] : &z, *v = i < b->n ? &b->a[i] : &z;
        LT_P_FREE(&r.a[i]);
        r.a[i] = LT_P_ADD(A->x, u, v, sign);
        LT_P_FREE(&z);
    }
    return LT_PFX(ap_take)(r.a, r.n);
}

static LT_PFX(AP) LT_PFX(ap_mul)(const LT_PFX(Alg)* A, const LT_PFX(AP)* a, const LT_PFX(AP)* b) {
    if (!a->n || !b->n) return LT_PFX(ap_zeros)(A, 0);
    LT_PFX(AP) r = LT_PFX(ap_zeros)(A, a->n + b->n - 1);
    for (int i = 0; i < a->n; i++)
        for (int j = 0; j < b->n; j++) {
            if (LT_P_DEG(&a->a[i]) < 0 || LT_P_DEG(&b->a[j]) < 0) continue;
            LT_P p = LT_PFX(a_mul)(A, &a->a[i], &b->a[j]), s = LT_P_ADD(A->x, &r.a[i + j], &p, 1);
            LT_P_FREE(&p); LT_P_FREE(&r.a[i + j]);
            r.a[i + j] = s;
        }
    return LT_PFX(ap_take)(r.a, r.n);
}

static LT_PFX(AP) LT_PFX(ap_one)(const LT_PFX(Alg)* A) {
    LT_P one = LT_P_ONE(A->x);
    LT_PFX(AP) p = LT_PFX(ap_term)(A, &one, 0);
    LT_P_FREE(&one);
    return p;
}

#define P LT_PFX(AP)
#define C LT_P
#define CTX const LT_PFX(Alg)*
#define PFX(name) LT_PFX(ra_##name)
#define KDIV LT_PFX(ADiv)
#define KBEZ LT_PFX(ABez)
#define KSF LT_PFX(ASf)
#define P_DEG(p) ((p)->n - 1)
#define P_LC(p) (&(p)->a[(p)->n - 1])
#define P_TERM(x, c, k) LT_PFX(ap_term)((x), (c), (k))
#define P_ADD(x, a, b, sign) LT_PFX(ap_add)((x), (a), (b), (sign))
#define P_MUL(x, a, b) LT_PFX(ap_mul)((x), (a), (b))
#define P_ONE(x) LT_PFX(ap_one)(x)
#define P_ZERO(x) LT_PFX(ap_zeros)((x), 0)
#define P_FREE(p) LT_PFX(ap_free)(p)
#define P_COEF(p, k) (&(p)->a[k])
#define P_LO(p) 0
#define P_N(p) ((p)->n)
#define P_TAKE(x, arr, lo, n) LT_PFX(ap_take)((arr), (n)) // lo is always 0 here
#define C_MUL(x, a, b) LT_PFX(a_mul)((x), (a), (b))
#define C_INV(x, a) LT_PFX(a_inv)((x), (a))
#define C_FREE(a) LT_P_FREE(a)
#define C_SCALE_INT(x, a, k) LT_PFX(a_scale_int)((x), (a), (k))
#define C_ADD(x, a, b, sign) LT_P_ADD((x)->x, (a), (b), (sign))
#define C_ZERO(x) LT_P_ZERO((x)->x)
#define C_IS_ZERO(a) (LT_P_DEG(a) < 0)
#define C_COPY(a) LT_P_COPY(NULL, (a)) // neither instance's copy reads its context
#include "kpoly_tmpl.h"
#undef P
#undef C
#undef CTX
#undef PFX
#undef KDIV
#undef KBEZ
#undef KSF
#undef P_DEG
#undef P_LC
#undef P_TERM
#undef P_ADD
#undef P_MUL
#undef P_ONE
#undef P_ZERO
#undef P_FREE
#undef P_COEF
#undef P_LO
#undef P_N
#undef P_TAKE
#undef C_MUL
#undef C_INV
#undef C_FREE
#undef C_SCALE_INT
#undef C_ADD
#undef C_ZERO
#undef C_IS_ZERO
#undef C_COPY

// --- between K(α)[θ] and K[θ] ---

// b ∈ K[θ] with constant coefficients in t
static LT_PFX(AP) LT_PFX(ap_embed)(const LT_PFX(Alg)* A, const LT_P* b) {
    int n = LT_P_DEG(b) + 1;
    LT_PFX(AP) p = LT_PFX(ap_zeros)(A, n);
    for (int k = 0; k < n; k++) {
        LT_C c = LT_P_COEFF(A->x, b, k);
        LT_P_FREE(&p.a[k]);
        p.a[k] = LT_P_TERM(A->x, &c, 0);
        LT_C_FREE(&c);
    }
    return LT_PFX(ap_take)(p.a, p.n);
}

// L_i ∈ K[θ] with g = Σ_i t^i·L_i: the t^i coefficients of g's θ-coefficients
static LT_P LT_PFX(ap_slice)(const LT_PFX(Alg)* A, const LT_PFX(AP)* g, int i) {
    LT_CTX x = A->x;
    (void)x;
    LT_P L = LT_P_ZERO(x);
    for (int k = 0; k < g->n; k++) {
        LT_C c = LT_P_COEFF(x, &g->a[k], i);
        LT_P t = LT_P_TERM(x, &c, k), s = LT_P_ADD(x, &L, &t, 1);
        LT_C_FREE(&c); LT_P_FREE(&t); LT_P_FREE(&L);
        L = s;
    }
    return L;
}

// Σ_i t^i·L[i] for n slices
static LT_PFX(AP) LT_PFX(ap_unslice)(const LT_PFX(Alg)* A, const LT_P* L) {
    LT_CTX x = A->x;
    (void)x;
    int m = 0;
    for (int i = 0; i < A->n; i++) if (LT_P_DEG(&L[i]) + 1 > m) m = LT_P_DEG(&L[i]) + 1;
    LT_PFX(AP) g = LT_PFX(ap_zeros)(A, m);
    for (int k = 0; k < m; k++)
        for (int i = 0; i < A->n; i++) {
            LT_C c = LT_P_COEFF(x, &L[i], k);
            LT_P t = LT_P_TERM(x, &c, i), s = LT_P_ADD(x, &g.a[k], &t, 1);
            LT_C_FREE(&c); LT_P_FREE(&t); LT_P_FREE(&g.a[k]);
            g.a[k] = s;
        }
    return LT_PFX(ap_take)(g.a, g.n);
}

// D(g) with α constant: D on each slice
static LT_PFX(AP) LT_PFX(ap_deriv)(const LT_PFX(Alg)* A, const LT_PFX(AP)* g) {
    LT_P* L = malloc((size_t)A->n * sizeof(LT_P));
    for (int i = 0; i < A->n; i++) {
        LT_P s = LT_PFX(ap_slice)(A, g, i);
        L[i] = LT_P_DERIV(A->x, &s);
        LT_P_FREE(&s);
    }
    LT_PFX(AP) r = LT_PFX(ap_unslice)(A, L);
    for (int i = 0; i < A->n; i++) LT_P_FREE(&L[i]);
    free(L);
    return r;
}

// Tr(g) = Σ_i ps[i]·L_i over the conjugates of α, ps[i] = Σ α^i
static LT_P LT_PFX(ap_trace)(const LT_PFX(Alg)* A, const LT_PFX(AP)* g) {
    LT_CTX x = A->x;
    (void)x;
    LT_P T = LT_P_ZERO(x);
    for (int i = 0; i < A->n; i++) {
        LT_P s = LT_PFX(ap_slice)(A, g, i);
        LT_C c = LT_C_RAT(x, &A->ps[i]);
        LT_P cs = LT_P_SCALE(x, &s, &c), u = LT_P_ADD(x, &T, &cs, 1);
        LT_P_FREE(&s); LT_P_FREE(&cs); LT_P_FREE(&T); LT_C_FREE(&c);
        T = u;
    }
    return T;
}

// --- f's data over Q ---

// f/lc(f): c[j] the coefficient of t^j, j <= n (c[n] = 1)
static Rational* LT_PFX(monic_coeffs)(const Polynomial* f, int n) {
    Rational* c = malloc((size_t)(n + 1) * sizeof(Rational)), lc;
    rat_init(&lc);
    poly_coeff(&lc, f, n);
    for (int j = 0; j <= n; j++) {
        rat_init(&c[j]);
        poly_coeff(&c[j], f, j);
        rat_div(&c[j], &c[j], &lc);
    }
    rat_clear(&lc);
    return c;
}

// the power sums ps[i] = Σ α^i, i < n, of monic c (Newton: p_i = -(i·c_{n-i} + Σ_{j<i} c_{n-j}·p_{i-j}))
static Rational* LT_PFX(power_sums)(const Rational* c, int n) {
    Rational* p = malloc((size_t)n * sizeof(Rational)), t;
    rat_init(&t);
    for (int i = 0; i < n; i++) {
        rat_init(&p[i]);
        if (i == 0) { mpq_set_si(p[0].val, n, 1); continue; }
        mpq_set_si(p[i].val, i, 1);
        rat_mul(&p[i], &p[i], &c[n - i]);
        for (int j = 1; j < i; j++) { rat_mul(&t, &c[n - j], &p[i - j]); rat_add(&p[i], &p[i], &t); }
        rat_neg(&p[i], &p[i]);
    }
    rat_clear(&t);
    return p;
}

static void LT_PFX(rats_free)(Rational* r, int n) {
    for (int i = 0; i < n; i++) rat_clear(&r[i]);
    free(r);
}

// Σ_j c[j]·t^j over K
static LT_P LT_PFX(poly_t)(LT_CTX x, const Rational* c, int n) {
    (void)x;
    LT_P F = LT_P_ZERO(x);
    for (int j = 0; j <= n; j++) {
        LT_C cj = LT_C_RAT(x, &c[j]);
        LT_P t = LT_P_TERM(x, &cj, j), s = LT_P_ADD(x, &F, &t, 1);
        LT_C_FREE(&cj); LT_P_FREE(&t); LT_P_FREE(&F);
        F = s;
    }
    return F;
}

// (a·b) mod d over K[θ]
static LT_P LT_PFX(mul_mod)(LT_CTX x, const LT_P* a, const LT_P* b, const LT_P* d) {
    (void)x;
    LT_P p = LT_P_MUL(x, a, b);
    LT_KDIV r = LT_P_DIVMOD(x, &p, d);
    LT_P_FREE(&p);
    LT_P_FREE(&r.q);
    return r.r;
}

// e = gcd(d, Σ_j c_j·h^j·D(d)^(n-j) mod d): homogeneous Horner, acc <- acc·h + c_j·D(d)^(n-j)
static LT_P LT_PFX(root_factor)(LT_CTX x, const LT_PFX(HD)* hd, const Rational* c, int n) {
    LT_P Dd = LT_P_DERIV(x, hd->d), acc = LT_P_ONE(x);
    LT_KDIV q = LT_P_DIVMOD(x, &Dd, hd->d);
    LT_P Dp = LT_P_COPY(x, &q.r);                          // D(d)^(n-j) mod d
    for (int j = n - 1; j >= 0; j--) {
        LT_C cj = LT_C_RAT(x, &c[j]);
        LT_P ah = LT_PFX(mul_mod)(x, &acc, hd->h, hd->d), cD = LT_P_SCALE(x, &Dp, &cj), s = LT_P_ADD(x, &ah, &cD, 1);
        LT_P Dn = LT_PFX(mul_mod)(x, &Dp, &q.r, hd->d);
        LT_C_FREE(&cj); LT_P_FREE(&ah); LT_P_FREE(&cD); LT_P_FREE(&acc); LT_P_FREE(&Dp);
        acc = s; Dp = Dn;
    }
    LT_KDIV r = LT_P_DIVMOD(x, &acc, hd->d);
    LT_P e = LT_P_GCD(x, hd->d, &r.r);
    LT_P_FREE(&Dd); LT_P_FREE(&acc); LT_P_FREE(&Dp);
    LT_KDIV_FREE(&q); LT_KDIV_FREE(&r);
    return e;
}

// --- the RootSum term ---

typedef struct {
    LT_CTX x;
    const Polynomial* f;    // irreducible over Q, degree n >= 3
    const LT_PFX(HD)* hd;   // h/d, d square-free and normal
} LT_PFX(RootIn);

typedef struct {
    LT_P* L;      // g_α = Σ_{i<n} t^i·L[i], L[i] ∈ K[θ] (monic in θ)
    int n;
    LT_FR deriv;  // D(Σ_α α·ln g_α), exactly
} LT_PFX(RootLog);

static void LT_PFX(root_log_free)(LT_PFX(RootLog)* r) {
    for (int i = 0; i < r->n; i++) LT_P_FREE(&r->L[i]);
    free(r->L);
    LT_FR_FREE(&r->deriv);
    r->L = NULL; r->n = 0;
}

// α = t, as a constant of K(α)[θ]
static LT_PFX(AP) LT_PFX(ap_alpha)(const LT_PFX(Alg)* A) {
    Rational one = rat_from_int(1);
    LT_C c = LT_C_RAT(A->x, &one);
    LT_P t = LT_P_TERM(A->x, &c, 1);
    LT_PFX(AP) a = LT_PFX(ap_term)(A, &t, 0);
    LT_P_FREE(&t); LT_C_FREE(&c);
    rat_clear(&one);
    return a;
}

// Tr(α·D(g)·(e/g)) / e: the derivative of Σ_α α·ln g_α
static LT_FR LT_PFX(root_deriv)(const LT_PFX(Alg)* A, const LT_PFX(AP)* g, const LT_P* e) {
    LT_PFX(AP) E = LT_PFX(ap_embed)(A, e), Dg = LT_PFX(ap_deriv)(A, g), cof = LT_PFX(ra_exact_div)(A, &E, g);
    LT_PFX(AP) a = LT_PFX(ap_alpha)(A), aDg = LT_PFX(ap_mul)(A, &a, &Dg), prod = LT_PFX(ap_mul)(A, &aDg, &cof);
    LT_FR d = LT_FR_MAKE(A->x, LT_PFX(ap_trace)(A, &prod), LT_P_COPY(A->x, e));
    LT_PFX(ap_free)(&E); LT_PFX(ap_free)(&Dg); LT_PFX(ap_free)(&cof);
    LT_PFX(ap_free)(&a); LT_PFX(ap_free)(&aDg); LT_PFX(ap_free)(&prod);
    return d;
}

// g_α = gcd(E, h - α·D(d) mod E) over K(α)[θ], monic
static LT_PFX(AP) LT_PFX(root_gcd)(const LT_PFX(Alg)* A, const LT_PFX(HD)* hd, const LT_P* e) {
    LT_P Dd = LT_P_DERIV(A->x, hd->d);
    LT_PFX(AP) E = LT_PFX(ap_embed)(A, e), H = LT_PFX(ap_embed)(A, hd->h), D = LT_PFX(ap_embed)(A, &Dd);
    LT_PFX(AP) a = LT_PFX(ap_alpha)(A), aD = LT_PFX(ap_mul)(A, &a, &D), W = LT_PFX(ap_add)(A, &H, &aD, -1);
    LT_PFX(ADiv) w = LT_PFX(ra_divmod)(A, &W, &E);
    LT_PFX(AP) g = LT_PFX(ra_gcd)(A, &E, &w.r);
    LT_P_FREE(&Dd);
    LT_PFX(ap_free)(&E); LT_PFX(ap_free)(&H); LT_PFX(ap_free)(&D); LT_PFX(ap_free)(&a);
    LT_PFX(ap_free)(&aD); LT_PFX(ap_free)(&W); LT_PFX(ap_free)(&w.q); LT_PFX(ap_free)(&w.r);
    return g;
}

// the term for f; false when a degree check fails (an internal error, not a proof)
static bool LT_PFX(root_log)(const LT_PFX(RootIn)* in, LT_PFX(RootLog)* out) {
    LT_CTX x = in->x;
    int n = polynomial_degree(in->f);
    Rational *c = LT_PFX(monic_coeffs)(in->f, n), *ps = LT_PFX(power_sums)(c, n);
    LT_PFX(Alg) A = { x, LT_PFX(poly_t)(x, c, n), n, ps };
    LT_P e = LT_PFX(root_factor)(x, in->hd, c, n);
    int de = LT_P_DEG(&e);
    bool ok = de > 0 && de % n == 0;                       // each of the n conjugates divides e
    *out = (LT_PFX(RootLog)){ NULL, 0, LT_FR_ZERO(x) };
    if (ok) {
        LT_PFX(AP) g = LT_PFX(root_gcd)(&A, in->hd, &e);
        ok = g.n - 1 == de / n;
        if (ok) {
            LT_FR_FREE(&out->deriv);
            out->deriv = LT_PFX(root_deriv)(&A, &g, &e);
            out->n = n;
            out->L = malloc((size_t)n * sizeof(LT_P));
            for (int i = 0; i < n; i++) out->L[i] = LT_PFX(ap_slice)(&A, &g, i);
        }
        LT_PFX(ap_free)(&g);
    }
    LT_P_FREE(&e); LT_P_FREE(&A.F);
    LT_PFX(rats_free)(c, n + 1); LT_PFX(rats_free)(ps, n);
    return ok;
}
