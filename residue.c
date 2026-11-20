// The residue test and the logarithmic part over K(θ) (Bronstein, Symbolic Integration I,
// §5.6 and §2.8). The residues of h/d are the roots of res_θ(d, h - z·D(d)), up to a factor in K
// the characteristic polynomial of s = h·D(d)⁻¹ mod d; they are constants iff every power sum
// Tr(s^k), k = 1 .. deg d, lies in Q (Newton's identities, exact). A nonconstant one proves
// ∫ h/d non-elementary (Thm 5.6.1). Otherwise R(z) ∈ Q[z] is factored: a linear factor gives
// c·ln(gcd(d, h - c·D(d))); a quadratic one gives the conjugate pair in real form, with the gcd
// taken over K(√M) (ksurd.c); degree >= 3 is not supported (it would need a RootSum over θ).
#include "residue.h"
#include "ksurd.h"
#include "factor.h"
#include <stdlib.h>

// --- the residue polynomial ---

// (a·b) mod d
static LPoly mul_mod(const LPoly* a, const LPoly* b, const LPoly* d) {
    LPoly p = lpoly_mul(a, b);
    KDiv q = kp_divmod(&p, d);
    lpoly_free(&p); lpoly_free(&q.q);
    return q.r;
}

// P[j] = Σ θ_i^j over the roots of monic d, j < deg d (Newton's identities)
static QFn* root_power_sums(const LPoly* d) {
    int n = kp_deg(d);
    QFn* P = malloc((size_t)n * sizeof(QFn));
    P[0] = qfn_const(n, 1);
    for (int j = 1; j < n; j++) {                         // P_j = -(j·c_{n-j} + Σ_{i<j} c_{n-i}·P_{j-i})
        QFn cj = kp_coeff(d, n - j), jq = qfn_const(j, 1), acc = qfn_mul(&jq, &cj);
        qfn_free(&jq);
        for (int i = 1; i < j; i++) {
            QFn ci = kp_coeff(d, n - i), t = qfn_mul(&ci, &P[j - i]), s = qfn_add(&acc, &t, 1);
            qfn_free(&acc); qfn_free(&ci); qfn_free(&t);
            acc = s;
        }
        QFn z = qfn_zero();
        P[j] = qfn_add(&z, &acc, -1);
        qfn_free(&z); qfn_free(&acc); qfn_free(&cj);
    }
    return P;
}

// Tr(a) on K[θ]/(d) for deg a < deg d: Σ a_j·P_j
static QFn trace(const LPoly* a, const QFn* P, int n) {
    QFn acc = qfn_zero();
    for (int j = 0; j < n; j++) {
        QFn aj = kp_coeff(a, j), t = qfn_mul(&aj, &P[j]), s = qfn_add(&acc, &t, 1);
        qfn_free(&aj); qfn_free(&t); qfn_free(&acc);
        acc = s;
    }
    return acc;
}

// the power sums p[k-1] = Tr(s^k), k = 1 .. deg d, of the residues; false when one is not
// a constant (the proof of non-elementarity)
static bool residue_power_sums(const KFrac* R, const Tower* t, Rational* p) {
    int n = kp_deg(&R->d);
    LPoly Dd = lpoly_derivative(&R->d, t), one = kp_one();
    KBezout z = kp_bezout(&Dd, &R->d, &one);              // z.s = D(d)⁻¹ mod d
    LPoly s = mul_mod(&R->n, &z.s, &R->d), pw = kp_one();
    QFn* P = root_power_sums(&R->d);
    bool ok = true;
    for (int k = 1; ok && k <= n; k++) {
        LPoly next = mul_mod(&pw, &s, &R->d);
        lpoly_free(&pw);
        pw = next;
        QFn tr = trace(&pw, P, n);
        ok = qfn_constant(&tr, &p[k - 1]);
        qfn_free(&tr);
    }
    for (int j = 0; j < n; j++) qfn_free(&P[j]);
    free(P);
    lpoly_free(&Dd); lpoly_free(&one); kbezout_free(&z); lpoly_free(&s); lpoly_free(&pw);
    return ok;
}

// the monic polynomial in z with power sums p[0..n-1]: k·e_k = Σ_{i=1}^{k} (-1)^(i-1)·e_{k-i}·p_i
static Polynomial* from_power_sums(const Rational* p, int n) {
    Rational* e = malloc((size_t)(n + 1) * sizeof(Rational));
    e[0] = rat_one_val();
    Polynomial* R = create_polynomial();
    add_term(R, &e[0], n);
    for (int k = 1; k <= n; k++) {
        rat_init(&e[k]);
        Rational t;
        rat_init(&t);
        for (int i = 1; i <= k; i++) {
            rat_mul(&t, &e[k - i], &p[i - 1]);
            if (i % 2) rat_add(&e[k], &e[k], &t); else rat_sub(&e[k], &e[k], &t);
        }
        Rational kk = rat_from_int(k);
        rat_div(&e[k], &e[k], &kk);
        rat_set(&t, &e[k]);                               // coefficient (-1)^k·e_k of z^(n-k)
        if (k % 2) rat_neg(&t, &t);
        add_term(R, &t, n - k);
        rat_clear(&kk); rat_clear(&t);
    }
    for (int k = 0; k <= n; k++) rat_clear(&e[k]);
    free(e);
    return R;
}

// --- collecting terms: rational residues c > 0, then c < 0, then the quadratic pairs ---

typedef struct {
    Expr* e[3];
    KFrac deriv;
    const Tower* t;
} Acc;

// adds e (taken) with derivative d (taken) to slot k
typedef struct {
    int slot;
    Expr* e;
    KFrac d;
} LogTerm;

static void acc_add(Acc* A, LogTerm* x) {
    A->e[x->slot] = expr_append_sum(A->e[x->slot], x->e);
    KFrac s = kfrac_add(&A->deriv, &x->d, A->t->kind);
    kfrac_free(&A->deriv); kfrac_free(&x->d);
    A->deriv = s;
}

// c·D(G)/G
static KFrac log_deriv(const Rational* c, const LPoly* G, const Tower* t) {
    QFn cq = qfn_rat(c);
    LPoly DG = lpoly_derivative(G, t);
    KFrac r = kfrac_make(kp_scale(&DG, &cq), kp_copy(G), t->kind);
    qfn_free(&cq); lpoly_free(&DG);
    return r;
}

// |c|·sqrt(m)·x with the sign outside (1 and sqrt(1) omitted; x may be NULL for a constant)
static Expr* surd_times(const Rational* c, const mpz_t m, Expr* x) {
    Rational a;
    rat_init(&a);
    rat_abs(&a, c);
    Expr* k = rat_is_one(&a) ? NULL : create_const_rat(&a);
    if (mpz_cmp_ui(m, 1) != 0) {
        Rational mr;
        rat_init(&mr);
        mpq_set_z(mr.val, m);
        Expr* s = create_sqrt(create_const_rat(&mr));
        k = k ? create_mul(k, s) : s;
        rat_clear(&mr);
    }
    Expr* e = !k ? (x ? x : create_const(1)) : x ? create_mul(k, x) : k;
    bool neg = rat_is_negative(c);
    rat_clear(&a);
    return neg ? create_neg(e) : e;
}

// the sign of T's leading coefficient (of the top θ power, then of the top x power)
static int lead_sign(const LPoly* T) {
    const Polynomial* top = T->a[T->n - 1].num;
    Rational c;
    rat_init(&c);
    poly_coeff(&c, top, polynomial_degree(top));
    int s = rat_is_negative(&c) ? -1 : 1;
    rat_clear(&c);
    return s;
}

// c·sqrt(m)·T with the coefficient of a one-term T = k·θ^j (k ∈ Q) folded in
static Expr* surd_poly(const Rational* c, const mpz_t m, const LPoly* T, const Tower* t) {
    Rational k;
    rat_init(&k);
    if (T->n != 1 || !qfn_constant(&T->a[0], &k)) { rat_clear(&k); return surd_times(c, m, lpoly_expr_order(T, t, true)); }
    rat_mul(&k, &k, c);
    QFn one = qfn_const(1, 1);
    Expr* e = surd_times(&k, m, T->lo ? lpoly_term_expr(&one, t, T->lo) : NULL);
    qfn_free(&one);
    rat_clear(&k);
    return e;
}

// -T
static LPoly kp_neg(const LPoly* T) {
    LPoly z = { 0, 0, NULL };
    return lpoly_add(&z, T, -1);
}

// c·ln(G)
static Expr* ln_expr(const Rational* c, const LPoly* G, const Tower* t) {
    mpz_t one;
    mpz_init_set_ui(one, 1);
    Expr* e = surd_times(c, one, create_ln(lpoly_expr_order(G, t, true)));
    mpz_clear(one);
    return e;
}

typedef struct {
    const KFrac* hd;
    const Tower* t;
} Ctx;

// c·ln(G_c) for G_c a primitive multiple of gcd(d, h - c·D(d))
static void rational_term(Acc* A, const Rational* c, const Ctx* X) {
    QFn cq = qfn_rat(c);
    LPoly Dd = lpoly_derivative(&X->hd->d, X->t), cDd = kp_scale(&Dd, &cq), w = lpoly_add(&X->hd->n, &cDd, -1);
    LPoly g = kp_gcd(&X->hd->d, &w), G = kp_primitive(&g);
    LogTerm x = { rat_is_negative(c) ? 1 : 0, ln_expr(c, &G, X->t), log_deriv(c, &G, X->t) };
    acc_add(A, &x);
    qfn_free(&cq); lpoly_free(&Dd); lpoly_free(&cDd); lpoly_free(&w); lpoly_free(&g); lpoly_free(&G);
}

// --- a conjugate pair α = a ± b·√M (M = m, or -m for complex α) ---

typedef struct {
    QuadRoots z;   // a = z.p, b = z.r, m = z.m
    Rational M;
    LPoly P, Q;    // g_α = P + √M·Q, scaled to coefficients in Q[x]
    const Tower* t;
} Pair;

// g_α = gcd(d, h - α·D(d)) over K(√M), with the x-denominators cleared
static void pair_gcd(Pair* q, const Ctx* X) {
    QFn a = qfn_rat(&q->z.p), b = qfn_rat(&q->z.r);
    LPoly Dd = lpoly_derivative(&X->hd->d, X->t), aDd = kp_scale(&Dd, &a), bDd = kp_scale(&Dd, &b), z = { 0, 0, NULL };
    LPoly wp = lpoly_add(&X->hd->n, &aDd, -1), wq = lpoly_add(&z, &bDd, -1);
    SPoly w = { wp, wq }, d = sp_from(&X->hd->d, &z);
    SPoly g = sp_gcd(&d, &w, &q->M);
    KPair c = kp_common_multiple(&g.p, &g.q);
    q->P = c.a; q->Q = c.b;
    sp_free(&w); sp_free(&d); sp_free(&g);
    qfn_free(&a); qfn_free(&b);
    lpoly_free(&Dd); lpoly_free(&aDd); lpoly_free(&bDd);
}

// P² - M·Q² (in K[θ])
static LPoly pair_norm(const Pair* q) {
    LPoly P2 = lpoly_mul(&q->P, &q->P), Q2 = lpoly_mul(&q->Q, &q->Q);
    QFn M = qfn_rat(&q->M);
    LPoly MQ2 = kp_scale(&Q2, &M), N = lpoly_add(&P2, &MQ2, -1);
    lpoly_free(&P2); lpoly_free(&Q2); lpoly_free(&MQ2); qfn_free(&M);
    return N;
}

// k·(P·D(Q) - Q·D(P))/(P² - M·Q²): the derivative of (k/2)/√M·ln|(P + √M·Q)/(P - √M·Q)|
static KFrac pair_deriv(const Pair* q, const Rational* k) {
    LPoly DP = lpoly_derivative(&q->P, q->t), DQ = lpoly_derivative(&q->Q, q->t);
    LPoly l = lpoly_mul(&q->P, &DQ), r = lpoly_mul(&q->Q, &DP), W = lpoly_add(&l, &r, -1);
    QFn kq = qfn_rat(k);
    KFrac f = kfrac_make(kp_scale(&W, &kq), pair_norm(q), q->t->kind);
    lpoly_free(&DP); lpoly_free(&DQ); lpoly_free(&l); lpoly_free(&r); lpoly_free(&W); qfn_free(&kq);
    return f;
}

// P + sign·sqrt(m)·Q as an expression
static Expr* p_pm_sqrt_q(const Pair* q, const LPoly* Q, int sign) {
    Rational s = rat_from_int(sign);
    Expr* e = create_add(lpoly_expr_order(&q->P, q->t, true), surd_poly(&s, q->z.m, Q, q->t));
    rat_clear(&s);
    return e;
}

// real α: b·√m·ln|P + √m·Q| - b·√m·ln|P - √m·Q|, odd in Q: written with lc(Q) > 0
static Expr* real_pair_expr(const Pair* q) {
    int sg = lead_sign(&q->Q);
    LPoly Q = sg > 0 ? kp_copy(&q->Q) : kp_neg(&q->Q);
    Rational b, nb;
    rat_init(&b); rat_init(&nb);
    rat_set(&b, &q->z.r);
    if (sg < 0) rat_neg(&b, &b);
    rat_neg(&nb, &b);
    Expr* e = surd_times(&b, q->z.m, create_ln(p_pm_sqrt_q(q, &Q, 1)));
    e = expr_append_sum(e, surd_times(&nb, q->z.m, create_ln(p_pm_sqrt_q(q, &Q, -1))));
    rat_clear(&b); rat_clear(&nb);
    lpoly_free(&Q);
    return e;
}

// --- Rioboo's LogToAtan with a scalar: 2·atan(κ·A/B), κ = ρ·√m, as Σ 2·atan(ρ_i·√m·T_i), T_i ∈ K[θ]
// (arguments polynomial in θ, so the sum has no jumps where B vanishes). With B·D - A·C' = G:
//   B | A:           2·atan(κ·A/B)
//   deg A < deg B:   2·atan(κ·A/B) = const - 2·atan(κ'·(-B)/A), κ' = 1/κ = √m/(ρ·m)
//   otherwise:       2·atan(κ·(A·D + B·C'/κ²)/G) + 2·atan(κ·D/C')     (Bronstein's LogToAtan with C = C'/κ)

typedef struct {
    Rational* rho;
    LPoly* T;
    int n;
} Atans;

static void atans_push(Atans* s, const Rational* rho, LPoly T) {
    s->rho = realloc(s->rho, (size_t)(s->n + 1) * sizeof(Rational));
    s->T = realloc(s->T, (size_t)(s->n + 1) * sizeof(LPoly));
    rat_init_set(&s->rho[s->n], rho);
    s->T[s->n++] = T;
}

static void atans_free(Atans* s) {
    for (int i = 0; i < s->n; i++) { rat_clear(&s->rho[i]); lpoly_free(&s->T[i]); }
    free(s->rho); free(s->T);
}

// the state 2·atan(ρ·√m·A/B)
typedef struct {
    LPoly A, B;
    Rational rho;
    const Rational* m;
} AtanState;

// one LogToAtan step; false once finished
static bool atan_step(AtanState* S, Atans* out) {
    KDiv qr = kp_divmod(&S->A, &S->B);
    bool done = qr.r.n == 0;
    if (done) { atans_push(out, &S->rho, qr.q); lpoly_free(&qr.r); return false; }
    kdiv_free(&qr);
    if (kp_deg(&S->A) < kp_deg(&S->B)) {                  // swap: (A, B, ρ) <- (-B, A, 1/(ρ·m))
        LPoly z = { 0, 0, NULL }, nb = lpoly_add(&z, &S->B, -1);
        lpoly_free(&S->B);
        S->B = S->A; S->A = nb;
        rat_mul(&S->rho, &S->rho, S->m);
        mpq_inv(S->rho.val, S->rho.val);
        return true;
    }
    LPoly G = kp_gcd(&S->A, &S->B), A1 = kp_exact_div(&S->A, &G), B1 = kp_exact_div(&S->B, &G), one = kp_one();
    KBezout z = kp_bezout(&B1, &A1, &one);                // D·B1 + t·A1 = 1: D = s, C' = -t
    LPoly zero = { 0, 0, NULL }, C = lpoly_add(&zero, &z.t, -1);
    Rational k2;                                          // 1/κ² = 1/(ρ²·m)
    rat_init(&k2);
    rat_mul(&k2, &S->rho, &S->rho);
    rat_mul(&k2, &k2, S->m);
    mpq_inv(k2.val, k2.val);
    QFn kq = qfn_rat(&k2);
    LPoly AD = lpoly_mul(&S->A, &z.s), BC = lpoly_mul(&S->B, &C), BCk = kp_scale(&BC, &kq), sum = lpoly_add(&AD, &BCk, 1);
    atans_push(out, &S->rho, kp_exact_div(&sum, &G));
    lpoly_free(&S->A); lpoly_free(&S->B);
    S->A = kp_copy(&z.s); S->B = C;
    lpoly_free(&G); lpoly_free(&A1); lpoly_free(&B1); lpoly_free(&one); kbezout_free(&z);
    lpoly_free(&AD); lpoly_free(&BC); lpoly_free(&BCk); lpoly_free(&sum);
    qfn_free(&kq); rat_clear(&k2);
    return true;
}

// Σ 2·ρ_i·D(T_i)/(1 + ρ_i²·m·T_i²) = 2·(P·D(Q) - Q·D(P))/(P² + m·Q²): the arctangents
// differentiate (divided by √m) to the pair's derivative, checked exactly
static bool atans_check(const Atans* s, const Pair* q) {
    TowerKind kind = q->t->kind;
    KFrac acc = { { 0, 0, NULL }, kp_one() };
    for (int i = 0; i < s->n; i++) {
        Rational k;
        rat_init(&k);
        rat_mul(&k, &s->rho[i], &s->rho[i]);
        rat_mul(&k, &k, &q->M);
        rat_neg(&k, &k);                                  // ρ²·m = -ρ²·M
        QFn kq = qfn_rat(&k), two = qfn_rat(&s->rho[i]);
        QFn tw = qfn_add(&two, &two, 1);
        LPoly T2 = lpoly_mul(&s->T[i], &s->T[i]), kT2 = kp_scale(&T2, &kq), one = kp_one(), den = lpoly_add(&one, &kT2, 1);
        LPoly DT = lpoly_derivative(&s->T[i], q->t);
        KFrac term = kfrac_make(kp_scale(&DT, &tw), den, kind), next = kfrac_add(&acc, &term, kind);
        kfrac_free(&acc); kfrac_free(&term);
        acc = next;
        lpoly_free(&T2); lpoly_free(&kT2); lpoly_free(&one); lpoly_free(&DT);
        qfn_free(&kq); qfn_free(&two); qfn_free(&tw); rat_clear(&k);
    }
    Rational two = rat_from_int(2);
    KFrac want = pair_deriv(q, &two), diff = kfrac_sub(&acc, &want, kind);
    bool same = diff.n.n == 0;
    kfrac_free(&acc); kfrac_free(&want); kfrac_free(&diff);
    rat_clear(&two);
    return same;
}

// complex α: -b·√m·Σ 2·atan(ρ_i·√m·T_i) from LogToAtan(√m·Q, P); NULL if the check fails
static Expr* complex_pair_expr(const Pair* q) {
    Rational m;
    rat_init(&m);
    mpq_set_z(m.val, q->z.m);
    AtanState S = { kp_copy(&q->Q), kp_copy(&q->P), rat_one_val(), &m };
    Atans s = { NULL, NULL, 0 };
    while (atan_step(&S, &s)) {}
    Expr* e = NULL;
    if (atans_check(&s, q)) {
        Rational c;                                       // every term: -2b·√m·atan(ρ·√m·T)
        rat_init(&c);
        rat_add(&c, &q->z.r, &q->z.r);
        rat_neg(&c, &c);
        for (int i = 0; i < s.n; i++) {                   // atan is odd: lc(ρ·T) > 0 inside
            int sg = lead_sign(&s.T[i]) * (rat_is_negative(&s.rho[i]) ? -1 : 1);
            LPoly T = sg > 0 ? kp_copy(&s.T[i]) : kp_neg(&s.T[i]);
            Rational ci;
            rat_init_set(&ci, &c);
            if (sg < 0) rat_neg(&ci, &ci);
            Expr* arg = surd_poly(&s.rho[i], q->z.m, &T, q->t);
            e = expr_append_sum(e, surd_times(&ci, q->z.m, create_atan(arg)));
            rat_clear(&ci); lpoly_free(&T);
        }
        rat_clear(&c);
    }
    lpoly_free(&S.A); lpoly_free(&S.B); rat_clear(&S.rho); rat_clear(&m);
    atans_free(&s);
    return e;
}

// a·ln|P² - M·Q²| plus the real or complex pair terms, for the roots of the quadratic factor f
static bool quadratic_term(Acc* A, const Polynomial* f, const Ctx* X) {
    Pair q;
    q.z = poly_quad_roots(f);
    q.t = X->t;
    rat_init(&q.M);
    mpq_set_z(q.M.val, q.z.m);
    if (q.z.is_complex) rat_neg(&q.M, &q.M);
    pair_gcd(&q, X);
    Expr* pe = q.z.is_complex ? complex_pair_expr(&q) : real_pair_expr(&q);
    if (pe) {
        Rational k;                                       // D(α ln g_α + ᾱ ln g_ᾱ) - a·D(N)/N = 2b·M·(P·DQ - Q·DP)/N
        rat_init(&k);
        rat_add(&k, &q.z.r, &q.z.r);
        rat_mul(&k, &k, &q.M);
        LogTerm x = { 2, pe, pair_deriv(&q, &k) };
        if (!rat_is_zero(&q.z.p)) {                       // a·ln|N|, first
            LPoly N = pair_norm(&q), Np = kp_primitive(&N);
            LogTerm l = { 2, ln_expr(&q.z.p, &Np, X->t), log_deriv(&q.z.p, &Np, X->t) };
            acc_add(A, &l);
            lpoly_free(&N); lpoly_free(&Np);
        }
        acc_add(A, &x);
        rat_clear(&k);
    }
    lpoly_free(&q.P); lpoly_free(&q.Q); rat_clear(&q.M); quad_roots_clear(&q.z);
    return pe != NULL;
}

// --- the entry point ---

// the terms for every irreducible factor of R(z) ∈ Q[z]
static ResidueStatus terms_for(const Polynomial* Rz, const Ctx* X, LogPart* out) {
    Factorization* fz = factorize(Rz);
    Acc A = { { NULL, NULL, NULL }, { { 0, 0, NULL }, kp_one() }, X->t };
    ResidueStatus st = RESIDUE_OK;
    for (int i = 0; st == RESIDUE_OK && i < fz->count; i++) {
        const Polynomial* f = fz->factors[i];
        int deg = polynomial_degree(f);
        if (deg == 1) {
            Rational a, b;                                // a·z + b: c = -b/a
            rat_init(&a); rat_init(&b);
            poly_coeff(&a, f, 1); poly_coeff(&b, f, 0);
            rat_div(&b, &b, &a);
            rat_neg(&b, &b);
            rational_term(&A, &b, X);
            rat_clear(&a); rat_clear(&b);
        } else if (deg > 2) {
            st = RESIDUE_UNSUPPORTED;
            out->reason = "residues that are algebraic numbers of degree 3 or more (a RootSum over exp, ln or tan) are not supported yet";
        } else if (!quadratic_term(&A, f, X)) {
            st = RESIDUE_UNSUPPORTED;
            out->reason = "internal check failed: the arctangents over exp/ln do not differentiate back";
        }
    }
    free_factorization(fz);
    out->e = NULL;
    for (int k = 0; k < 3; k++) {
        if (st == RESIDUE_OK && A.e[k]) out->e = expr_append_sum(out->e, A.e[k]);
        else free_expr(A.e[k]);
    }
    out->deriv = A.deriv;
    return st;
}

ResidueStatus residue_log_part(const KFrac* hd, const Tower* t, LogPart* out) {
    *out = (LogPart){ NULL, { { 0, 0, NULL }, kp_one() }, NULL };
    int n = hd->n.n ? kp_deg(&hd->d) : 0;                 // no residues when h = 0
    if (n == 0) return RESIDUE_OK;
    Rational* p = malloc((size_t)n * sizeof(Rational));
    for (int k = 0; k < n; k++) rat_init(&p[k]);
    ResidueStatus st = RESIDUE_NONCONSTANT;
    if (residue_power_sums(hd, t, p)) {
        Polynomial* Rz = from_power_sums(p, n);
        kfrac_free(&out->deriv);
        st = terms_for(Rz, &(Ctx){ hd, t }, out);
        free_polynomial(Rz);
    }
    for (int k = 0; k < n; k++) rat_clear(&p[k]);
    free(p);
    return st;
}

void log_part_free(LogPart* L) {
    free_expr(L->e);
    kfrac_free(&L->deriv);
    L->e = NULL;
}
