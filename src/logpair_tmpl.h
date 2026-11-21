// The algebra of a conjugate pair of quadratic residues α = a ± b·√M over a field K and θ
// (residue.c over K = Q(x), risch2.c over K1 = Q(x)(θ1)), written once and instantiated per
// type at compile time like kpoly_tmpl.h. Include this file after defining
//   LT_P, LT_C, LT_FR      polynomials in θ, their coefficients in K, and reduced fractions
//   LT_CTX                 the derivation context passed to every operation
//   LT_PFX(name)           the prefix for the generated names
//   LT_KDIV {q, r}, LT_KBEZ {s, t}, with LT_KDIV_FREE(d), LT_KBEZ_FREE(z)
//   coefficients: LT_C_ZERO(x), LT_C_RAT(x, r) (r ∈ Q), LT_C_ADD(x, a, b, sign), LT_C_MUL(x, a, b),
//                 LT_C_INV(x, a), LT_C_FREE(a)
//   polynomials:  LT_P_ZERO(x), LT_P_ONE(x), LT_P_COPY(x, a), LT_P_FREE(a), LT_P_DEG(a) (-1 for 0),
//                 LT_P_COEFF(x, a, k) (a copy, 0 outside), LT_P_TERM(x, c, k) (c·θ^k, c copied),
//                 LT_P_ADD(x, a, b, sign), LT_P_MUL(x, a, b), LT_P_SCALE(x, a, c) (c ∈ K),
//                 LT_P_DERIV(x, a) (the derivation D), LT_P_DIVMOD, LT_P_GCD (monic),
//                 LT_P_EXACT_DIV, LT_P_BEZOUT(x, a, b, c) (s·a + t·b = c)
//   fractions:    LT_FR_MAKE(x, n, d) (takes both), LT_FR_ADD(x, a, b, sign), LT_FR_ZERO(x),
//                 LT_FR_IS_ZERO(f), LT_FR_FREE(f)
// The generated functions: pair_make (the gcd g_α = gcd(d, h - α·D(d)) = P + √M·Q over K(√M)),
// pair_norm, pair_deriv, pair_free, and log_to_atan (Rioboo's arctangents, checked exactly).
// rootsum_tmpl.h follows with the same macros; logterm_undef.h then undefines them.

// --- K(√M)[θ]: just enough for a monic gcd ---

// u + √M·w
typedef struct {
    LT_C u, w;
} LT_PFX(SElt);

// p + √M·q
typedef struct {
    LT_P p, q;
} LT_PFX(SPoly);

typedef struct {
    LT_CTX x;
    const Rational* M;
} LT_PFX(SCtx);

static int LT_PFX(sp_deg)(const LT_PFX(SPoly)* a) {
    int dp = LT_P_DEG(&a->p), dq = LT_P_DEG(&a->q);
    return dp > dq ? dp : dq;
}

static void LT_PFX(sp_free)(LT_PFX(SPoly)* a) { LT_P_FREE(&a->p); LT_P_FREE(&a->q); }
static void LT_PFX(se_free)(LT_PFX(SElt)* e) { LT_C_FREE(&e->u); LT_C_FREE(&e->w); }

// (u1 + √M·w1)(u2 + √M·w2) = (u1·u2 + M·w1·w2) + √M·(u1·w2 + w1·u2)
static LT_PFX(SElt) LT_PFX(se_mul)(const LT_PFX(SCtx)* s, const LT_PFX(SElt)* a, const LT_PFX(SElt)* b) {
    LT_CTX x = s->x;
    (void)x;
    LT_C uu = LT_C_MUL(x, &a->u, &b->u), ww = LT_C_MUL(x, &a->w, &b->w), uw = LT_C_MUL(x, &a->u, &b->w), wu = LT_C_MUL(x, &a->w, &b->u);
    LT_C m = LT_C_RAT(x, s->M), mww = LT_C_MUL(x, &m, &ww);
    LT_PFX(SElt) r = { LT_C_ADD(x, &uu, &mww, 1), LT_C_ADD(x, &uw, &wu, 1) };
    LT_C_FREE(&uu); LT_C_FREE(&ww); LT_C_FREE(&uw); LT_C_FREE(&wu); LT_C_FREE(&m); LT_C_FREE(&mww);
    return r;
}

// 1/(u + √M·w) = (u - √M·w)/(u² - M·w²); the norm is nonzero since M is not a square
static LT_PFX(SElt) LT_PFX(se_inv)(const LT_PFX(SCtx)* s, const LT_PFX(SElt)* a) {
    LT_CTX x = s->x;
    (void)x;
    LT_C uu = LT_C_MUL(x, &a->u, &a->u), ww = LT_C_MUL(x, &a->w, &a->w), m = LT_C_RAT(x, s->M), mww = LT_C_MUL(x, &m, &ww);
    LT_C n = LT_C_ADD(x, &uu, &mww, -1), in = LT_C_INV(x, &n), z = LT_C_ZERO(x), nw = LT_C_ADD(x, &z, &a->w, -1);
    LT_PFX(SElt) r = { LT_C_MUL(x, &a->u, &in), LT_C_MUL(x, &nw, &in) };
    LT_C_FREE(&uu); LT_C_FREE(&ww); LT_C_FREE(&m); LT_C_FREE(&mww); LT_C_FREE(&n); LT_C_FREE(&in); LT_C_FREE(&z); LT_C_FREE(&nw);
    return r;
}

static LT_PFX(SElt) LT_PFX(sp_coeff)(LT_CTX x, const LT_PFX(SPoly)* a, int k) {
    (void)x;
    return (LT_PFX(SElt)){ LT_P_COEFF(x, &a->p, k), LT_P_COEFF(x, &a->q, k) };
}

// c·θ^k
typedef struct {
    LT_PFX(SElt) c;
    int k;
} LT_PFX(STerm);

// a·c·θ^k
static LT_PFX(SPoly) LT_PFX(sp_mul_term)(const LT_PFX(SCtx)* s, const LT_PFX(SPoly)* a, const LT_PFX(STerm)* t) {
    LT_CTX x = s->x;
    (void)x;
    LT_P cu = LT_P_TERM(x, &t->c.u, t->k), cw = LT_P_TERM(x, &t->c.w, t->k);
    LT_P pu = LT_P_MUL(x, &a->p, &cu), qw = LT_P_MUL(x, &a->q, &cw), pw = LT_P_MUL(x, &a->p, &cw), qu = LT_P_MUL(x, &a->q, &cu);
    LT_C m = LT_C_RAT(x, s->M);
    LT_P mqw = LT_P_SCALE(x, &qw, &m);
    LT_PFX(SPoly) r = { LT_P_ADD(x, &pu, &mqw, 1), LT_P_ADD(x, &pw, &qu, 1) };
    LT_P_FREE(&cu); LT_P_FREE(&cw); LT_P_FREE(&pu); LT_P_FREE(&qw); LT_P_FREE(&pw); LT_P_FREE(&qu); LT_P_FREE(&mqw);
    LT_C_FREE(&m);
    return r;
}

static LT_PFX(SPoly) LT_PFX(sp_copy)(LT_CTX x, const LT_PFX(SPoly)* a) {
    (void)x;
    return (LT_PFX(SPoly)){ LT_P_COPY(x, &a->p), LT_P_COPY(x, &a->q) };
}

// a mod b, from the top one term at a time
static LT_PFX(SPoly) LT_PFX(sp_rem)(const LT_PFX(SCtx)* s, const LT_PFX(SPoly)* a, const LT_PFX(SPoly)* b) {
    LT_CTX x = s->x;
    LT_PFX(SPoly) r = LT_PFX(sp_copy)(x, a);
    int db = LT_PFX(sp_deg)(b);
    LT_PFX(SElt) lb = LT_PFX(sp_coeff)(x, b, db), ib = LT_PFX(se_inv)(s, &lb);
    while (LT_PFX(sp_deg)(&r) >= db) {
        int dr = LT_PFX(sp_deg)(&r);
        LT_PFX(SElt) lr = LT_PFX(sp_coeff)(x, &r, dr);
        LT_PFX(STerm) t = { LT_PFX(se_mul)(s, &lr, &ib), dr - db };
        LT_PFX(SPoly) tb = LT_PFX(sp_mul_term)(s, b, &t);
        LT_PFX(SPoly) n = { LT_P_ADD(x, &r.p, &tb.p, -1), LT_P_ADD(x, &r.q, &tb.q, -1) };
        LT_PFX(sp_free)(&r); LT_PFX(sp_free)(&tb); LT_PFX(se_free)(&lr); LT_PFX(se_free)(&t.c);
        r = n;
    }
    LT_PFX(se_free)(&lb); LT_PFX(se_free)(&ib);
    return r;
}

// monic gcd over K(√M); a or b nonzero
static LT_PFX(SPoly) LT_PFX(sp_gcd)(const LT_PFX(SCtx)* s, const LT_PFX(SPoly)* a, const LT_PFX(SPoly)* b) {
    LT_CTX x = s->x;
    LT_PFX(SPoly) r0 = LT_PFX(sp_copy)(x, a), r1 = LT_PFX(sp_copy)(x, b);
    while (LT_PFX(sp_deg)(&r1) >= 0) {                    // (r0, r1) <- (r1, r0 mod r1)
        LT_PFX(SPoly) r = LT_PFX(sp_rem)(s, &r0, &r1);
        LT_PFX(sp_free)(&r0);
        r0 = r1; r1 = r;
    }
    LT_PFX(SElt) lc = LT_PFX(sp_coeff)(x, &r0, LT_PFX(sp_deg)(&r0));
    LT_PFX(STerm) t = { LT_PFX(se_inv)(s, &lc), 0 };
    LT_PFX(SPoly) g = LT_PFX(sp_mul_term)(s, &r0, &t);
    LT_PFX(se_free)(&lc); LT_PFX(se_free)(&t.c);
    LT_PFX(sp_free)(&r0); LT_PFX(sp_free)(&r1);
    return g;
}

// --- the pair ---

typedef struct {
    QuadRoots z;   // a = z.p, b = z.r, m = z.m
    Rational M;    // m, or -m for complex α
    LT_P P, Q;     // g_α = P + √M·Q
    LT_CTX x;
} LT_PFX(Pair);

// h/d with d square-free
typedef struct {
    const LT_P *h, *d;
} LT_PFX(HD);

// gcd(d, (u² - M·v²) mod d) for w = u + √M·v: the product g_α·g_ᾱ, found without surds, so the gcd
// over K(√M) runs against this small factor of d instead of d itself (Euclid there swells fast)
static LT_P LT_PFX(norm_factor)(LT_CTX x, const LT_PFX(SPoly)* w, const LT_P* d, const Rational* M) {
    (void)x;
    LT_P u2 = LT_P_MUL(x, &w->p, &w->p), v2 = LT_P_MUL(x, &w->q, &w->q);
    LT_C m = LT_C_RAT(x, M);
    LT_P mv2 = LT_P_SCALE(x, &v2, &m), n = LT_P_ADD(x, &u2, &mv2, -1);
    LT_KDIV r = LT_P_DIVMOD(x, &n, d);
    LT_P g = LT_P_GCD(x, d, &r.r);
    LT_P_FREE(&u2); LT_P_FREE(&v2); LT_P_FREE(&mv2); LT_P_FREE(&n); LT_C_FREE(&m);
    LT_KDIV_FREE(&r);
    return g;
}

// the roots of the quadratic f and g_α = gcd(d, h - α·D(d)), monic over K(√M)
static void LT_PFX(pair_make)(LT_PFX(Pair)* q, const Polynomial* f, const LT_PFX(HD)* hd) {
    LT_CTX x = q->x;
    q->z = poly_quad_roots(f);
    rat_init(&q->M);
    mpq_set_z(q->M.val, q->z.m);
    if (q->z.is_complex) rat_neg(&q->M, &q->M);
    LT_C a = LT_C_RAT(x, &q->z.p), b = LT_C_RAT(x, &q->z.r);
    LT_P Dd = LT_P_DERIV(x, hd->d), aDd = LT_P_SCALE(x, &Dd, &a), bDd = LT_P_SCALE(x, &Dd, &b), z = LT_P_ZERO(x);
    LT_PFX(SPoly) w = { LT_P_ADD(x, hd->h, &aDd, -1), LT_P_ADD(x, &z, &bDd, -1) };
    LT_PFX(SPoly) d = { LT_PFX(norm_factor)(x, &w, hd->d, &q->M), LT_P_ZERO(x) };
    LT_PFX(SPoly) g = LT_PFX(sp_gcd)(&(LT_PFX(SCtx)){ x, &q->M }, &d, &w);
    q->P = g.p; q->Q = g.q;
    LT_PFX(sp_free)(&w); LT_PFX(sp_free)(&d);
    LT_C_FREE(&a); LT_C_FREE(&b);
    LT_P_FREE(&Dd); LT_P_FREE(&aDd); LT_P_FREE(&bDd); LT_P_FREE(&z);
}

static void LT_PFX(pair_free)(LT_PFX(Pair)* q) {
    LT_P_FREE(&q->P); LT_P_FREE(&q->Q);
    rat_clear(&q->M);
    quad_roots_clear(&q->z);
}

// P² - M·Q² (in K[θ])
static LT_P LT_PFX(pair_norm)(const LT_PFX(Pair)* q) {
    LT_CTX x = q->x;
    (void)x;
    LT_P P2 = LT_P_MUL(x, &q->P, &q->P), Q2 = LT_P_MUL(x, &q->Q, &q->Q);
    LT_C M = LT_C_RAT(x, &q->M);
    LT_P MQ2 = LT_P_SCALE(x, &Q2, &M), N = LT_P_ADD(x, &P2, &MQ2, -1);
    LT_P_FREE(&P2); LT_P_FREE(&Q2); LT_P_FREE(&MQ2); LT_C_FREE(&M);
    return N;
}

// k·(P·D(Q) - Q·D(P))/(P² - M·Q²): the derivative of (k/2)/√M·ln|(P + √M·Q)/(P - √M·Q)|
static LT_FR LT_PFX(pair_deriv)(const LT_PFX(Pair)* q, const Rational* k) {
    LT_CTX x = q->x;
    LT_P DP = LT_P_DERIV(x, &q->P), DQ = LT_P_DERIV(x, &q->Q);
    LT_P l = LT_P_MUL(x, &q->P, &DQ), r = LT_P_MUL(x, &q->Q, &DP), W = LT_P_ADD(x, &l, &r, -1);
    LT_C kq = LT_C_RAT(x, k);
    LT_FR f = LT_FR_MAKE(x, LT_P_SCALE(x, &W, &kq), LT_PFX(pair_norm)(q));
    LT_P_FREE(&DP); LT_P_FREE(&DQ); LT_P_FREE(&l); LT_P_FREE(&r); LT_P_FREE(&W); LT_C_FREE(&kq);
    return f;
}

// --- Rioboo's LogToAtan with a scalar: 2·atan(κ·A/B), κ = ρ·√m, as Σ 2·atan(ρ_i·√m·T_i), T_i ∈ K[θ]
// (arguments polynomial in θ, so the sum has no jumps where B vanishes). With B·D - A·C' = G:
//   B | A:           2·atan(κ·A/B)
//   deg A < deg B:   2·atan(κ·A/B) = const - 2·atan(κ'·(-B)/A), κ' = 1/κ = √m/(ρ·m)
//   otherwise:       2·atan(κ·(A·D + B·C'/κ²)/G) + 2·atan(κ·D/C')     (Bronstein's LogToAtan with C = C'/κ)

typedef struct {
    Rational* rho;
    LT_P* T;
    int n;
} LT_PFX(Atans);

static void LT_PFX(atans_push)(LT_PFX(Atans)* s, const Rational* rho, LT_P T) {
    s->rho = realloc(s->rho, (size_t)(s->n + 1) * sizeof(Rational));
    s->T = realloc(s->T, (size_t)(s->n + 1) * sizeof(LT_P));
    rat_init_set(&s->rho[s->n], rho);
    s->T[s->n++] = T;
}

static void LT_PFX(atans_free)(LT_PFX(Atans)* s) {
    for (int i = 0; i < s->n; i++) { rat_clear(&s->rho[i]); LT_P_FREE(&s->T[i]); }
    free(s->rho); free(s->T);
}

// the state 2·atan(ρ·√m·A/B)
typedef struct {
    LT_P A, B;
    Rational rho;
    const Rational* m;
    LT_CTX x;
} LT_PFX(AtanState);

// -a
static LT_P LT_PFX(p_neg)(LT_CTX x, const LT_P* a) {
    (void)x;
    LT_P z = LT_P_ZERO(x), r = LT_P_ADD(x, &z, a, -1);
    LT_P_FREE(&z);
    return r;
}

// one LogToAtan step; false once finished
static bool LT_PFX(atan_step)(LT_PFX(AtanState)* S, LT_PFX(Atans)* out) {
    LT_CTX x = S->x;
    LT_KDIV qr = LT_P_DIVMOD(x, &S->A, &S->B);
    if (LT_P_DEG(&qr.r) < 0) { LT_PFX(atans_push)(out, &S->rho, qr.q); LT_P_FREE(&qr.r); return false; }
    LT_KDIV_FREE(&qr);
    if (LT_P_DEG(&S->A) < LT_P_DEG(&S->B)) {              // swap: (A, B, ρ) <- (-B, A, 1/(ρ·m))
        LT_P nb = LT_PFX(p_neg)(x, &S->B);
        LT_P_FREE(&S->B);
        S->B = S->A; S->A = nb;
        rat_mul(&S->rho, &S->rho, S->m);
        mpq_inv(S->rho.val, S->rho.val);
        return true;
    }
    LT_P G = LT_P_GCD(x, &S->A, &S->B), A1 = LT_P_EXACT_DIV(x, &S->A, &G), B1 = LT_P_EXACT_DIV(x, &S->B, &G), one = LT_P_ONE(x);
    LT_KBEZ z = LT_P_BEZOUT(x, &B1, &A1, &one);           // D·B1 + t·A1 = 1: D = s, C' = -t
    LT_P C = LT_PFX(p_neg)(x, &z.t);
    Rational k2;                                          // 1/κ² = 1/(ρ²·m)
    rat_init(&k2);
    rat_mul(&k2, &S->rho, &S->rho);
    rat_mul(&k2, &k2, S->m);
    mpq_inv(k2.val, k2.val);
    LT_C kq = LT_C_RAT(x, &k2);
    LT_P AD = LT_P_MUL(x, &S->A, &z.s), BC = LT_P_MUL(x, &S->B, &C), BCk = LT_P_SCALE(x, &BC, &kq), sum = LT_P_ADD(x, &AD, &BCk, 1);
    LT_PFX(atans_push)(out, &S->rho, LT_P_EXACT_DIV(x, &sum, &G));
    LT_P_FREE(&S->A); LT_P_FREE(&S->B);
    S->A = LT_P_COPY(x, &z.s); S->B = C;
    LT_P_FREE(&G); LT_P_FREE(&A1); LT_P_FREE(&B1); LT_P_FREE(&one); LT_KBEZ_FREE(&z);
    LT_P_FREE(&AD); LT_P_FREE(&BC); LT_P_FREE(&BCk); LT_P_FREE(&sum);
    LT_C_FREE(&kq); rat_clear(&k2);
    return true;
}

// 2ρ·D(T)/(1 + ρ²·m·T²), m = -M: the derivative of 2·atan(ρ·√m·T) divided by √m
static LT_FR LT_PFX(atan_deriv)(const LT_PFX(Pair)* q, const Rational* rho, const LT_P* T) {
    LT_CTX x = q->x;
    Rational k, r2;
    rat_init(&k); rat_init(&r2);
    rat_mul(&k, rho, rho);
    rat_mul(&k, &k, &q->M);
    rat_neg(&k, &k);                                      // ρ²·m = -ρ²·M
    rat_add(&r2, rho, rho);
    LT_C kq = LT_C_RAT(x, &k), tw = LT_C_RAT(x, &r2);
    LT_P T2 = LT_P_MUL(x, T, T), kT2 = LT_P_SCALE(x, &T2, &kq), one = LT_P_ONE(x), den = LT_P_ADD(x, &one, &kT2, 1);
    LT_P DT = LT_P_DERIV(x, T);
    LT_FR f = LT_FR_MAKE(x, LT_P_SCALE(x, &DT, &tw), den);
    LT_P_FREE(&T2); LT_P_FREE(&kT2); LT_P_FREE(&one); LT_P_FREE(&DT);
    LT_C_FREE(&kq); LT_C_FREE(&tw);
    rat_clear(&k); rat_clear(&r2);
    return f;
}

// Σ 2·ρ_i·D(T_i)/(1 + ρ_i²·m·T_i²) = 2·(P·D(Q) - Q·D(P))/(P² + m·Q²): the arctangents
// differentiate (divided by √m) to the pair's derivative, checked exactly
static bool LT_PFX(atans_check)(const LT_PFX(Atans)* s, const LT_PFX(Pair)* q) {
    LT_CTX x = q->x;
    (void)x;
    LT_FR acc = LT_FR_ZERO(x);
    for (int i = 0; i < s->n; i++) {
        LT_FR term = LT_PFX(atan_deriv)(q, &s->rho[i], &s->T[i]), next = LT_FR_ADD(x, &acc, &term, 1);
        LT_FR_FREE(&acc); LT_FR_FREE(&term);
        acc = next;
    }
    Rational two = rat_from_int(2);
    LT_FR want = LT_PFX(pair_deriv)(q, &two), diff = LT_FR_ADD(x, &acc, &want, -1);
    bool same = LT_FR_IS_ZERO(&diff);
    LT_FR_FREE(&acc); LT_FR_FREE(&want); LT_FR_FREE(&diff);
    rat_clear(&two);
    return same;
}

// complex α (M = -m < 0): LogToAtan(√m·Q, P) as Σ 2·atan(ρ_i·√m·T_i) in out; false if the check fails
static bool LT_PFX(log_to_atan)(const LT_PFX(Pair)* q, LT_PFX(Atans)* out) {
    LT_CTX x = q->x;
    Rational m;
    rat_init(&m);
    mpq_set_z(m.val, q->z.m);
    LT_PFX(AtanState) S = { LT_P_COPY(x, &q->Q), LT_P_COPY(x, &q->P), rat_one_val(), &m, x };
    *out = (LT_PFX(Atans)){ NULL, NULL, 0 };
    while (LT_PFX(atan_step)(&S, out)) {}
    bool ok = LT_PFX(atans_check)(out, q);
    LT_P_FREE(&S.A); LT_P_FREE(&S.B); rat_clear(&S.rho); rat_clear(&m);
    return ok;
}
