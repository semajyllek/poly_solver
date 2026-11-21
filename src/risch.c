// The Risch algorithm for one exponential: θ = exp(v), v ∈ Q[x] nonconstant, integrand
// f = Σ a_k·θ^k with a_k ∈ Q(x). By Liouville's theorem (and θ transcendental over Q(x)),
// ∫ f is elementary iff ∫ f = ∫ a_0 + Σ_{k≠0} y_k·θ^k with y_k ∈ Q(x), where each y_k solves
// the Risch differential equation y' + k·v'·y = a_k. Here k·v' is a polynomial, so:
//   * y has poles only at poles of a_k, and ν_p(y) = ν_p(a_k) + 1 there (y' dominates).
//     A simple pole of a_k therefore has no solution: ∫ exp(x)/x is not elementary.
//   * otherwise y = q/h with h = gcd(D, D') (D = den a_k) and q a polynomial solving
//     h·q' + (k·v'·h - h')·q = N·h²/D, whose leading term fixes deg q = deg C - deg B;
//     q follows term by term from the top, and a nonzero remainder means no solution.
// The same file handles θ = ln(u), denominators in θ (Hermite reduction and the residue test),
// and θ = tan(H) with its special factor 1 + θ² (see the sections below).
#include "risch.h"
#include "expfield.h"
#include "integrate.h"
#include "modgcd.h"
#include "kpoly.h"
#include "residue.h"
#include "expsum.h"
#include "tanfield.h"
#include "rischk1.h"
#include "risch2.h"
#include <stdlib.h>

// --- the Risch differential equation y' + f·y = g, f a nonzero polynomial ---

typedef struct {
    const Polynomial* f; // k·v'
    const QFn* g;        // a_k, in lowest terms
} Rde;

// q with A·q' + B·q = C, deg B > deg A - 1 (so the equation is triangular from the top);
// NULL when there is none
static Polynomial* solve_poly_rde(const Polynomial* A, const Polynomial* B, const Polynomial* C) {
    Polynomial *q = create_polynomial(), *R = poly_copy(C);
    Rational lb, lr;
    rat_init(&lb); rat_init(&lr);
    poly_coeff(&lb, B, polynomial_degree(B));
    while (!poly_is_zero(R) && polynomial_degree(R) >= polynomial_degree(B)) {
        poly_coeff(&lr, R, polynomial_degree(R));
        rat_div(&lr, &lr, &lb);
        Polynomial* t = create_polynomial();             // the next term of q
        add_term(t, &lr, polynomial_degree(R) - polynomial_degree(B));
        Polynomial *dt = poly_derivative(t), *At = poly_multiply(A, dt), *Bt = poly_multiply(B, t);
        Polynomial *s = poly_add(At, Bt), *nr = poly_subtract(R, s), *nq = poly_add(q, t);
        free_polynomial(dt); free_polynomial(At); free_polynomial(Bt); free_polynomial(s);
        free_polynomial(R); free_polynomial(q); free_polynomial(t);
        R = nr; q = nq;
    }
    bool solved = poly_is_zero(R);
    free_polynomial(R);
    rat_clear(&lb); rat_clear(&lr);
    if (!solved) { free_polynomial(q); return NULL; }
    return q;
}

// y ∈ Q(x) with y' + f·y = g, or false when there is none
static bool solve_rde(const Rde* r, QFn* y) {
    if (qfn_is_zero(r->g)) { *y = qfn_zero(); return true; }
    const Polynomial *D = r->g->den, *N = r->g->num;
    Polynomial* Dp = poly_derivative(D);
    Polynomial* h = poly_gcd_fast(D, Dp);                // Π p^(e-1) over D = Π p^e
    Polynomial* h2 = poly_multiply(h, h);
    PolyDivResult qd = poly_divmod(h2, D);                // D | h² iff every pole has order >= 2
    bool ok = poly_is_zero(qd.remainder);
    Polynomial* q = NULL;
    if (ok) {
        Polynomial *C = poly_multiply(N, qd.quotient), *hp = poly_derivative(h), *fh = poly_multiply(r->f, h);
        Polynomial* B = poly_subtract(fh, hp);
        q = solve_poly_rde(h, B, C);
        free_polynomial(C); free_polynomial(hp); free_polynomial(fh); free_polynomial(B);
    }
    if (q) *y = qfn_from(q, poly_copy(h));
    free_polynomial(Dp); free_polynomial(h2);
    free_polynomial(qd.quotient); free_polynomial(qd.remainder);
    free_polynomial(h);
    return q != NULL;
}

// --- integration ---

RischResult failure(RischStatus s, const char* why) { return (RischResult){ s, NULL, why }; }

// the constant 0, as a constant or as the zero polynomial
bool is_zero_expr(const Expr* e) {
    return (e->type == EXPR_CONST && rat_is_zero(&e->constant)) || (e->type == EXPR_POLY && poly_is_zero(e->poly));
}

// Σ_{k≠0} y_k·θ^k with y_k' + k·v'·y_k = a_k; false (with *bad the failing index) when some
// equation has no solution
bool exp_part(const ExpForm* F, LPoly* Y, int* bad) {
    Y->lo = F->f.lo;
    Y->n = F->f.n;
    Y->a = malloc((size_t)(Y->n > 0 ? Y->n : 1) * sizeof(QFn));
    Polynomial* vp = poly_derivative(F->t.v);
    bool ok = true;
    int i = 0;
    for (; ok && i < Y->n; i++) {
        int k = F->f.lo + i;
        if (k == 0) { Y->a[i] = qfn_zero(); continue; }
        Rational kk = rat_from_int(k);
        Polynomial* f = poly_scale(vp, &kk);
        ok = solve_rde(&(Rde){ f, &F->f.a[i] }, &Y->a[i]);
        free_polynomial(f);
        rat_clear(&kk);
    }
    if (!ok) { *bad = i - 1; Y->n = i - 1; lpoly_free(Y); } // the failed slot was never filled
    free_polynomial(vp);
    return ok;
}

// the θ-free part a_0 of f, as an expression (0 when absent)
static Expr* a0_expr(const LPoly* f) {
    int i = -f->lo;
    if (i < 0 || i >= f->n) return create_const(0);
    return qfn_expr(&f->a[i]);
}

// d/dx (Σ y_k θ^k) = Σ_{k≠0} a_k θ^k exactly (the y_k are unique, so this is a full check)
static bool check_exp_part(const ExpForm* F, const LPoly* Y) {
    LPoly dY = lpoly_derivative(Y, &F->t);
    QFn zero = qfn_zero();
    LPoly target = { F->f.lo, F->f.n, malloc((size_t)(F->f.n > 0 ? F->f.n : 1) * sizeof(QFn)) };
    for (int i = 0; i < F->f.n; i++) target.a[i] = F->f.lo + i == 0 ? qfn_copy(&zero) : qfn_copy(&F->f.a[i]);
    bool same = lpoly_equal(&dY, &target);
    lpoly_free(&dY); lpoly_free(&target); qfn_free(&zero);
    return same;
}

// --- the log case: θ = ln(u), L = u'/u ---
//
// For f = Σ_{k=0}^{n} a_k·θ^k, Liouville's theorem gives ∫ f = Σ_{k=1}^{n+1} b_k·θ^k + (an integral
// over Q(x)) with b_{n+1} constant and b_k ∈ Q(x). Matching θ^k in the derivative:
//   a_k = b_k' + (k+1)·b_{k+1}·L.
// Write b_k = r_k + d_k (d_k constant). From the top, g_k = a_k - (k+1)·r_{k+1}·L must satisfy
// ∫ g_k = r_k + (k+1)·d_{k+1}·θ: rational plus a constant multiple of ln u and nothing else
// ("limited integration"). If some level has other logs (or atans), ∫ f is not elementary.
// Level 0 is free: ∫ (a_0 - r_1·L) is any rational integral.

// ∫ g = R + c·θ with R ∈ Q(x), c ∈ Q: Hermite gives ∫ g = R + ∫ A/D, D square-free, and then
// A/D must equal c·L exactly; false otherwise (also when Hermite's self-check fails)
static bool limited_integrate(const QFn* g, const QFn* L, QFn* R, Rational* c) {
    HermiteSplit h;
    bool ok = hermite_split(g->num, g->den, &h);
    QFn rem = qfn_from(poly_copy(h.a), poly_copy(h.dstar));
    rat_set_si(c, 0, 1);
    if (ok && !qfn_is_zero(&rem)) {                       // rem = c·L: compare rem·L.den·... exactly
        QFn ratio = qfn_from(poly_multiply(rem.num, L->den), poly_multiply(rem.den, L->num));
        ok = polynomial_degree(ratio.num) <= 0 && polynomial_degree(ratio.den) == 0;
        if (ok) poly_coeff(c, ratio.num, 0);
        qfn_free(&ratio);
    }
    if (ok) *R = qfn_from(poly_copy(h.rnum), poly_copy(h.rden));
    qfn_free(&rem);
    hermite_split_free(&h);
    return ok;
}

// c·q for a rational constant c
static QFn qfn_scaled(const QFn* q, const Rational* c) {
    return qfn_from(poly_scale(q->num, c), poly_copy(q->den));
}

bool solve_log_levels(const ExpForm* F, LogSolution* S) {
    QFn L = tower_L(&F->t);
    int n = S->n;
    QFn r_next = qfn_zero();                              // r_{k+1}
    Rational d_next;                                      // d_{k+1}, set at level k
    rat_init(&d_next);
    bool ok = true;
    for (int k = n; ok && k >= 1; k--) {
        QFn a = kp_coeff(&F->f, k), rL = qfn_mul(&r_next, &L);
        Rational kk = rat_from_int(k + 1);
        QFn t = qfn_scaled(&rL, &kk), g = qfn_add(&a, &t, -1);
        QFn R;
        Rational c;
        rat_init(&c);
        ok = limited_integrate(&g, &L, &R, &c);
        if (ok) {
            rat_div(&d_next, &c, &kk);                    // (k+1)·d_{k+1} = c
            Rational one = rat_one_val();
            QFn dq = qfn_from(poly_constant(&d_next), poly_constant(&one));
            rat_clear(&one);
            qfn_free(&S->b[k + 1]);
            S->b[k + 1] = qfn_add(&r_next, &dq, 1);       // b_{k+1} = r_{k+1} + d_{k+1}
            qfn_free(&dq);
            qfn_free(&r_next);
            r_next = R;
        }
        qfn_free(&a); qfn_free(&rL); qfn_free(&t); qfn_free(&g);
        rat_clear(&kk); rat_clear(&c);
    }
    if (ok) {
        qfn_free(&S->b[1]);
        S->b[1] = qfn_copy(&r_next);                      // d_1 = 0: level 0 absorbs it
        QFn a0 = kp_coeff(&F->f, 0), rL = qfn_mul(&r_next, &L);
        qfn_free(&S->g0);
        S->g0 = qfn_add(&a0, &rL, -1);
        qfn_free(&a0); qfn_free(&rL);
    }
    qfn_free(&r_next);
    rat_clear(&d_next);
    qfn_free(&L);
    return ok;
}

// Σ_{k>=1} b_k θ^k as a polynomial in θ
LPoly log_part_poly(const LogSolution* S) {
    LPoly B = { 0, 0, NULL };
    for (int k = 1; k <= S->n + 1; k++) {
        LPoly t = lpoly_term(qfn_copy(&S->b[k]), k), sum = lpoly_add(&B, &t, 1);
        lpoly_free(&B); lpoly_free(&t);
        B = sum;
    }
    return B;
}

static RischResult integrate_log(const ExpForm* F, const Expr* f) {
    LogSolution S = { NULL, qfn_zero(), F->f.n > 0 ? F->f.lo + F->f.n - 1 : 0 };
    S.b = malloc((size_t)(S.n + 2) * sizeof(QFn));
    for (int k = 0; k <= S.n + 1; k++) S.b[k] = qfn_zero();
    RischResult r = { RISCH_OK, NULL, NULL };
    if (!solve_log_levels(F, &S))
        r = (RischResult){ RISCH_NONELEMENTARY, expr_copy(f), "has no elementary antiderivative (Risch: a ln level needs a logarithm other than the one present)" };
    LPoly B = log_part_poly(&S);
    if (r.status == RISCH_OK) {                           // exact check: B' + g0 = f
        LPoly dB = lpoly_derivative(&B, &F->t), G = lpoly_term(qfn_copy(&S.g0), 0), s = lpoly_add(&dB, &G, 1);
        if (!lpoly_equal(&s, &F->f)) r = failure(RISCH_UNSUPPORTED, "internal check failed: the ln part does not differentiate back");
        lpoly_free(&dB); lpoly_free(&G); lpoly_free(&s);
    }
    if (r.status == RISCH_OK) {
        Expr* g0 = qfn_expr(&S.g0);
        Expr* I0 = integrate(g0, "x");                    // the free level 0 (self-checked)
        free_expr(g0);
        if (!I0) r = failure(RISCH_UNSUPPORTED, integrate_error());
        else {
            Expr* eb = lpoly_expr(&B, &F->t);
            r.result = is_zero_expr(I0) ? (free_expr(I0), eb) : is_zero_expr(eb) ? (free_expr(eb), I0) : expr_append_sum(eb, I0);
        }
    }
    lpoly_free(&B);
    for (int k = 0; k <= S.n + 1; k++) qfn_free(&S.b[k]);
    free(S.b);
    qfn_free(&S.g0);
    return r;
}

// Σ_{k≠0} y_k·θ^k + ∫ a_0 for a Laurent polynomial f in θ = exp(v)
static RischResult integrate_exp_poly(const ExpForm* F) {
    LPoly Y;
    int bad = 0;
    if (!exp_part(F, &Y, &bad))
        return (RischResult){ RISCH_NONELEMENTARY, lpoly_term_expr(&F->f.a[bad], &F->t, F->f.lo + bad),
                              "has no elementary antiderivative (Risch: y' + k·v'·y = a has no rational solution)" };
    RischResult r = { RISCH_OK, NULL, NULL };
    if (!check_exp_part(F, &Y)) r = failure(RISCH_UNSUPPORTED, "internal check failed: the exp part does not differentiate back");
    Expr* a0 = a0_expr(&F->f);
    Expr* I0 = r.status == RISCH_OK ? integrate(a0, "x") : NULL;   // rational part (self-checked)
    if (r.status == RISCH_OK && !I0) r = failure(RISCH_UNSUPPORTED, integrate_error());
    if (r.status == RISCH_OK) {
        Expr* ey = lpoly_expr(&Y, &F->t);
        r.result = is_zero_expr(I0) ? (free_expr(I0), ey) : is_zero_expr(ey) ? (free_expr(ey), I0) : expr_append_sum(ey, I0);
    }
    else free_expr(I0);
    free_expr(a0);
    lpoly_free(&Y);
    return r;
}

// the integral of a polynomial (Laurent for exp) in θ
RischResult integrate_poly(const ExpForm* F, const Expr* f) {
    return F->t.kind == TOWER_LOG ? integrate_log(F, f) : integrate_exp_poly(F);
}

// --- denominators in θ: Hermite reduction and the residue test ---
//
// For f = N/E (E monic in θ, and θ ∤ E for exp, where θ | θ' makes θ special):
//   1. Bezout splits f = p + b/E with p a (Laurent) polynomial in θ and deg b < deg E.
//   2. Hermite reduction over K[θ], K = Q(x): square-free parts V of E by d/dθ, each step
//      solving B·U·D(V) + C·V = -A/j with the derivation D (every factor of E is normal, so
//      gcd(V, D(V)) = 1), gives ∫ b/E = g + ∫ h/d with d square-free and deg h < deg d.
//   3. Rothstein–Trager: the residues of h/d are the roots of res_θ(d, h - z·D(d)), up to a
//      factor in K the characteristic polynomial of s = h·D(d)⁻¹ mod d. ∫ f is elementary only
//      if these roots are constants (Bronstein, Symbolic Integration I, Thm 5.6.1), i.e. the
//      power sums Tr(s^k), k = 1 .. deg d, all lie in Q; then ∫ h/d = Σ c·ln(g_c) + (an element
//      of K[θ]), g_c = gcd(d, h - c·D(d)) over the distinct roots c.
//   4. f - D(g) - Σ c·D(g_c)/g_c is then a polynomial in θ (checked exactly), integrated as above.


// b with f = p + b/den, deg b < deg den: for f = N·θ^-m (N a polynomial), Bezout gives
// s·θ^m + t·den = N, so f = s/den + t·θ^-m
LPoly proper_part(const ExpForm* F) {
    int m = F->f.n && F->f.lo < 0 ? -F->f.lo : 0;
    LPoly N = kp_shift(&F->f, m), one = kp_one(), tm = kp_shift(&one, m);
    KBezout z = kp_bezout(&tm, &F->den, &N);
    LPoly b = z.s;
    lpoly_free(&z.t); lpoly_free(&N); lpoly_free(&one); lpoly_free(&tm);
    return b;
}

// (n/d)·a
static LPoly kp_times_const(const LPoly* a, long n, long d) {
    QFn c = qfn_const(n, d);
    LPoly r = kp_scale(a, &c);
    qfn_free(&c);
    return r;
}

typedef struct {
    const LPoly *U, *V, *UDV; // U = den/V^i and U·D(V)
    int j;
} HermiteStep;

// A/(U·V^(j+1)) = D(B/V^j) + (-j·C - U·D(B))/(U·V^j) with B·U·D(V) + C·V = -A/j
static void hermite_step(Reduced* R, LPoly* A, const HermiteStep* st, const Tower* t) {
    LPoly rhs = kp_times_const(A, -1, st->j);
    KBezout z = kp_bezout(st->UDV, st->V, &rhs);
    KFrac term = kfrac_make(kp_copy(&z.s), kp_pow(st->V, st->j), t->kind), g = kfrac_add(&R->g, &term, t->kind);
    kfrac_free(&R->g);
    R->g = g;
    LPoly DB = lpoly_derivative(&z.s, t), UDB = lpoly_mul(st->U, &DB), jC = kp_times_const(&z.t, -st->j, 1);
    lpoly_free(A);
    *A = lpoly_add(&jC, &UDB, -1);
    lpoly_free(&rhs); kbezout_free(&z); kfrac_free(&term);
    lpoly_free(&DB); lpoly_free(&UDB); lpoly_free(&jC);
}

// ∫ b/den = g + ∫ h/d (Bronstein's quadratic HermiteReduce with the derivation D)
Reduced hermite_reduce(const LPoly* b, const ExpForm* F) {
    Reduced R = { { { 0, 0, NULL }, kp_one() }, { 0, 0, NULL }, { 0, 0, NULL } };
    KSquareFree sf = kp_square_free(&F->den);
    LPoly A = kp_copy(b), den = kp_copy(&F->den);
    for (int i = 2; i <= sf.m; i++) {
        const LPoly* V = &sf.f[i - 1];
        if (kp_deg(V) <= 0) continue;
        LPoly Vi = kp_pow(V, i), U = kp_exact_div(&den, &Vi);
        LPoly DV = lpoly_derivative(V, &F->t), UDV = lpoly_mul(&U, &DV);
        for (int j = i - 1; j >= 1; j--) hermite_step(&R, &A, &(HermiteStep){ &U, V, &UDV, j }, &F->t);
        lpoly_free(&den);
        den = lpoly_mul(&U, V);                           // V^i -> V
        lpoly_free(&Vi); lpoly_free(&U); lpoly_free(&DV); lpoly_free(&UDV);
    }
    KDiv q = kp_divmod(&A, &den);                         // in the exp case A may have grown
    R.h = q.r;
    R.d = den;
    lpoly_free(&q.q); lpoly_free(&A);
    ksf_free(&sf);
    return R;
}

void reduced_free(Reduced* R) { kfrac_free(&R->g); lpoly_free(&R->h); lpoly_free(&R->d); }

// the integrand f = F (in its tower F->t) and its reduction R
typedef struct {
    const ExpForm* F;
    const Reduced* R;
    const Expr* f;
} FracCtx;

// f - D(g) - D(log part)
static KFrac remainder_frac(const ExpForm* F, const Reduced* R, const LogPart* L) {
    TowerKind kind = F->t.kind;
    KFrac f = { kp_copy(&F->f), kp_copy(&F->den) }, Dg = kfrac_derivative(&R->g, &F->t);
    KFrac a = kfrac_sub(&f, &Dg, kind), acc = kfrac_sub(&a, &L->deriv, kind);
    kfrac_free(&f); kfrac_free(&Dg); kfrac_free(&a);
    return acc;
}

// the remainder as a polynomial in θ, or false when it is not one (an internal failure: by
// Rothstein–Trager it always is)
static bool remainder_poly(const FracCtx* X, const LogPart* L, LPoly* r) {
    KFrac acc = remainder_frac(X->F, X->R, L);
    bool ok = kp_is_one(&acc.d);
    if (ok) *r = kp_copy(&acc.n);
    kfrac_free(&acc);
    return ok;
}

// acc + e, acc possibly zero; takes both
Expr* plus(Expr* acc, Expr* e) {
    if (is_zero_expr(e)) { free_expr(e); return acc; }
    if (is_zero_expr(acc)) { free_expr(acc); return e; }
    return expr_append_sum(acc, e);
}

// the polynomial remainder's integral, then g, then the log part (taken)
static RischResult finish_fraction(const FracCtx* X, LogPart* L) {
    const ExpForm* F = X->F;
    LPoly r;
    if (!remainder_poly(X, L, &r)) return failure(RISCH_UNSUPPORTED, "internal check failed: the fraction in exp/ln does not reduce to a polynomial");
    ExpForm P = { F->t, r, kp_one(), NULL };              // shares F's tower
    RischResult rest = integrate_poly(&P, X->f);
    lpoly_free(&P.f); lpoly_free(&P.den);
    if (rest.status != RISCH_OK) return rest;
    Expr* acc = rest.result;                              // x - ln(exp(x) + 1): the polynomial part first
    if (X->R->g.n.n) acc = plus(acc, kfrac_expr(&X->R->g, &F->t));
    if (L->e) { acc = plus(acc, L->e); L->e = NULL; }
    rest.result = acc;
    return rest;
}

#define RESIDUES_NOT_CONSTANT "has no elementary antiderivative (Risch: its residues are not constants)"

RischResult integrate_fraction(const ExpForm* F, const Expr* f) {
    LPoly b = proper_part(F);
    Reduced R = hermite_reduce(&b, F);
    KFrac hd = { R.h, R.d };
    LogPart L;
    RischResult r = { RISCH_OK, NULL, NULL };
    switch (residue_log_part(&hd, &F->t, &L)) {
        case RESIDUE_NONCONSTANT: r = (RischResult){ RISCH_NONELEMENTARY, kfrac_expr(&hd, &F->t), RESIDUES_NOT_CONSTANT }; break;
        case RESIDUE_UNSUPPORTED: r = failure(RISCH_UNSUPPORTED, L.reason); break;
        case RESIDUE_OK: r = finish_fraction(&(FracCtx){ F, &R, f }, &L); break;
    }
    log_part_free(&L);
    reduced_free(&R);
    lpoly_free(&b);
    return r;
}

// --- θ = tan(H): the special factor 1 + θ² (tanfield.h) ---
//
// 1 + θ² divides its own derivative (D(1 + θ²) = 2H'·θ·(1 + θ²)), so for f = N/D with
// D = (1 + θ²)^m·E, gcd(E, 1 + θ²) = 1 (Bronstein, Symbolic Integration I, §5.10):
//   1. f = P + S/(1 + θ²)^m + A/E with P a polynomial, deg S < 2m and deg A < deg E.
//   2. A/E is normal: Hermite reduction and the residue test as above.
//   3. r = f - D(g) - D(logs) lies in K[θ, 1/(1 + θ²)] (checked exactly). Its polynomial part
//      drops to p1·θ + p0 by D(c·θ^(n-1)) = c'·θ^(n-1) + (n-1)·c·H'·(θ^(n-2) + θ^n), and
//      p1·θ = c·D(ln(1 + θ²)) with c = p1/(2H'), which must be a constant (in z = e^{2iH},
//      θ = -i(z - 1)/(z + 1): c is the residue at z = -1, up to a factor -2).
//   4. The rest, p0 + S/(1 + θ²)^m with θ^a/(1 + θ²)^m = sin(H)^a·cos(H)^(2m-a), is a Laurent
//      polynomial in z, integrated by expsum.c: one RDE per power of z, ∫ p0 by the rational
//      integrator. Parts 2-4 are independent: D maps proper fractions in z (z ∤ den) to
//      proper fractions, so only the special part has nonzero powers of z.

// D = (1 + θ²)^m·E with gcd(E, 1 + θ²) = 1
typedef struct {
    LPoly S, E;   // S = (1 + θ²)^m
    int m;
} SpecialSplit;

static SpecialSplit split_special(const LPoly* D) {
    SpecialSplit z = { kp_one(), kp_copy(D), 0 };
    LPoly q = tan_special();
    for (;;) {
        KDiv d = kp_divmod(&z.E, &q);
        bool divides = d.r.n == 0;
        if (divides) { lpoly_free(&z.E); z.E = kp_copy(&d.q); LPoly s = lpoly_mul(&z.S, &q); lpoly_free(&z.S); z.S = s; z.m++; }
        kdiv_free(&d);
        if (!divides) break;
    }
    lpoly_free(&q);
    return z;
}

// H' as an element of K
static QFn tower_dv(const Tower* t) {
    QFn h = qfn_poly(t->v), d = qfn_derivative(&h);
    qfn_free(&h);
    return d;
}

// q with P - D(q) = p1·θ + p0 (P is replaced by p1·θ + p0): the top term p_n·θ^n is
// D(c·θ^(n-1))'s, c = p_n/((n-1)·H')
static LPoly reduce_tan_poly(LPoly* P, const Tower* t) {
    LPoly q = { 0, 0, NULL };
    QFn d = tower_dv(t);
    for (int n = kp_deg(P); n >= 2; n = kp_deg(P)) {
        QFn k = qfn_const(n - 1, 1), kh = qfn_mul(&k, &d), ik = qfn_inv(&kh), lc = kp_coeff(P, n);
        LPoly term = lpoly_term(qfn_mul(&lc, &ik), n - 1), Dt = lpoly_derivative(&term, t);
        LPoly nP = lpoly_add(P, &Dt, -1), nq = lpoly_add(&q, &term, 1);
        lpoly_free(P); lpoly_free(&q);
        *P = nP; q = nq;
        lpoly_free(&term); lpoly_free(&Dt);
        qfn_free(&k); qfn_free(&kh); qfn_free(&ik); qfn_free(&lc);
    }
    qfn_free(&d);
    return q;
}

// sin(H)^a·cos(H)^b (either power may be 0)
static Expr* sin_cos_power(const Tower* t, int a, int b) {
    Expr* e = NULL;
    const int pw[2] = { a, b };
    for (int i = 0; i < 2; i++) {
        if (pw[i] == 0) continue;
        Expr* f = create_poly_expr(poly_copy(t->v));
        f = i == 0 ? create_sin(f) : create_cos(f);
        if (pw[i] > 1) f = create_pow(f, create_const(pw[i]));
        e = e ? create_mul(e, f) : f;
    }
    return e ? e : create_const(1);
}

// Σ s_a·sin(H)^a·cos(H)^(n-a) for S = Σ s_a·θ^a, i.e. cos(H)^n·S(tan H); NULL for 0
static Expr* homogeneous_expr(const LPoly* S, int n, const Tower* t) {
    Expr* acc = NULL;
    for (int i = 0; i < S->n; i++) {
        if (qfn_is_zero(&S->a[i])) continue;
        int a = S->lo + i;
        Expr* term = create_mul(qfn_expr(&S->a[i]), sin_cos_power(t, a, n - a));
        acc = acc ? create_add(acc, term) : term;
    }
    return acc;
}

// the remainder r = N/(1 + θ²)^m as P (polynomial) and S/(1 + θ²)^m; false when its
// denominator is not a power of 1 + θ² (an internal failure)
typedef struct {
    LPoly P, S;
    int m;
} TanRest;

static bool tan_rest(const KFrac* r, TanRest* out) {
    SpecialSplit z = split_special(&r->d);
    bool ok = kp_is_one(&z.E);
    if (ok) {
        KDiv d = kp_divmod(&r->n, &z.S);
        *out = (TanRest){ d.q, d.r, z.m };
    }
    lpoly_free(&z.S); lpoly_free(&z.E);
    return ok;
}

// p0 + S/(1 + θ²)^m as an expression in sin(H) and cos(H) (NULL for 0)
static Expr* special_expr(const QFn* p0, const TanRest* R, const Tower* t) {
    Expr *a = qfn_is_zero(p0) ? NULL : qfn_expr(p0), *b = homogeneous_expr(&R->S, 2 * R->m, t);
    return !a ? b : !b ? a : create_add(a, b);
}

#define TAN_P1_NOT_CONSTANT "has no elementary antiderivative (Risch: p·tan(H) needs p/(2H') constant)"

// c = p1/(2H') when it is a rational constant
static bool log_special_coeff(const QFn* p1, const Tower* t, Rational* c) {
    QFn d = tower_dv(t), two = qfn_const(2, 1), d2 = qfn_mul(&two, &d), id2 = qfn_inv(&d2), cq = qfn_mul(p1, &id2);
    bool ok = qfn_constant(&cq, c);
    qfn_free(&d); qfn_free(&two); qfn_free(&d2); qfn_free(&id2); qfn_free(&cq);
    return ok;
}

// c·ln(1 + θ²) = -2c·ln|cos(H)|
static Expr* log_cos_term(const Rational* c, const Tower* t) {
    Rational k = rat_from_int(-2);
    rat_mul(&k, &k, c);
    Expr* e = create_mul(create_const_rat(&k), create_ln(create_cos(create_poly_expr(poly_copy(t->v)))));
    rat_clear(&k);
    return e;
}

// ∫ (P + S/(1 + θ²)^m), the part after the normal one: ∫ special + q + c·ln(1 + θ²), added to *acc
static RischResult integrate_tan_rest(TanRest* R, const Tower* t, Expr** acc) {
    LPoly q = reduce_tan_poly(&R->P, t);
    QFn p1 = kp_coeff(&R->P, 1), p0 = kp_coeff(&R->P, 0);
    Rational c;
    rat_init(&c);
    RischResult r = { RISCH_OK, NULL, NULL };
    if (!log_special_coeff(&p1, t, &c)) r = (RischResult){ RISCH_NONELEMENTARY, lpoly_term_expr(&p1, t, 1), TAN_P1_NOT_CONSTANT };
    else {
        Expr* sp = special_expr(&p0, R, t);
        if (sp) r = integrate_expsum(sp);                 // self-checked
        free_expr(sp);
        if (r.status == RISCH_OK) {
            if (r.result) *acc = plus(*acc, r.result);
            *acc = plus(*acc, lpoly_expr(&q, t));
            if (mpq_sgn(c.val)) *acc = plus(*acc, log_cos_term(&c, t));
            r.result = NULL;
        }
    }
    lpoly_free(&q);
    qfn_free(&p1); qfn_free(&p0);
    rat_clear(&c);
    return r;
}

// the normal part s/E of f = N/D (D = (1 + θ²)^m·E): N mod D = s·(1 + θ²)^m + t·E, deg s < deg E
static LPoly normal_numerator(const ExpForm* F, const SpecialSplit* z) {
    KDiv qd = kp_divmod(&F->f, &F->den);
    KBezout bz = kp_bezout(&z->S, &z->E, &qd.r);
    LPoly s = kp_copy(&bz.s);
    kdiv_free(&qd); kbezout_free(&bz);
    return s;
}

// --- the answer's logs in sin and cos: ln|G(tan H)| = ln|cos(H)^n·G(tan H)| - n·ln|cos(H)| ---

// one term e of the answer; when e = c·ln(arg), arg (else NULL), and G is arg as a nonconstant
// polynomial in θ when it is one (else G.n = 0)
typedef struct {
    const Expr* e;
    const Expr* arg;
    Rational c;
    LPoly G;
} LogTerm;

typedef struct {
    LogTerm* t;
    int n;
} LogTerms;

// e = c·ln(arg) for ln, -ln, c·ln, -(c·ln)
static const Expr* match_log(const Expr* e, Rational* c) {
    if (e->type == EXPR_NEG) { const Expr* a = match_log(e->operand, c); if (a) rat_neg(c, c); return a; }
    if (e->type == EXPR_LN) { rat_set_si(c, 1, 1); return e->ln_arg; }
    if (e->type == EXPR_BINOP && e->binop.op == OP_MUL && e->binop.left->type == EXPR_CONST && e->binop.right->type == EXPR_LN) {
        rat_set(c, &e->binop.left->constant);
        return e->binop.right->ln_arg;
    }
    return NULL;
}

// the terms of a sum (answers are built with + only; any other node is one term)
static void collect_terms(const Expr* e, LogTerms* T) {
    if (e->type == EXPR_BINOP && e->binop.op == OP_ADD) {
        collect_terms(e->binop.left, T);
        collect_terms(e->binop.right, T);
        return;
    }
    T->t = realloc(T->t, (size_t)(T->n + 1) * sizeof(LogTerm));
    LogTerm* L = &T->t[T->n++];
    L->e = e;
    L->G = (LPoly){ 0, 0, NULL };
    rat_init(&L->c);
    L->arg = match_log(e, &L->c);
}

static void log_terms_free(LogTerms* T) {
    for (int i = 0; i < T->n; i++) { rat_clear(&T->t[i].c); lpoly_free(&T->t[i].G); }
    free(T->t);
}

// cos(H) itself
static bool is_cos_h(const Expr* arg, const Tower* t) {
    if (arg->type != EXPR_COS) return false;
    Expr* c = canonicalize((Expr*)arg->fn_arg);
    bool same = c && c->type == EXPR_POLY;
    if (same) { Polynomial* d = poly_subtract(c->poly, t->v); same = poly_is_zero(d); free_polynomial(d); }
    free_expr(c);
    return same;
}

// c·ln|Ĝ| with Ĝ = cos(H)^n·G(tan H) in canonical form, squared when n is odd and H = B/2 (so
// that it is written with B); *out = NULL when Ĝ is a constant. False when Ĝ has no canonical
// form (the caller then keeps the tan form)
static bool rewritten_log(const LogTerm* L, const Tower* t, Expr** out) {
    int n = kp_deg(&L->G);
    bool square = t->half && n % 2;
    Expr* g = homogeneous_expr(&L->G, n, t);
    if (square) g = create_pow(g, create_const(2));
    const char* why = NULL;
    Expr* s = simplify_expsum(g, &why);
    free_expr(g);
    *out = NULL;
    if (!s) return false;
    if (s->type == EXPR_CONST || (s->type == EXPR_POLY && polynomial_degree(s->poly) <= 0)) { free_expr(s); return true; }
    char* str = expr_to_string(s);
    if (str[0] == '-') { Expr* m = create_neg(s); s = simplify_expsum(m, &why); free_expr(m); } // ln|-g| = ln|g|
    free(str);
    if (!s) return false;
    Rational k;
    rat_init(&k);
    rat_set(&k, &L->c);
    if (square) { Rational h = rat_from_int(2); rat_div(&k, &k, &h); rat_clear(&h); }
    *out = create_mul(create_const_rat(&k), create_ln(s));
    rat_clear(&k);
    return true;
}

// the ln|cos(H)| terms cancel once every ln|G(tan H)| is rewritten
static bool logs_cancel(const LogTerms* T, const Tower* t) {
    Rational k, m;
    rat_init(&k); rat_init(&m);
    for (int i = 0; i < T->n; i++) {
        const LogTerm* L = &T->t[i];
        if (!L->arg) continue;
        if (is_cos_h(L->arg, t)) rat_add(&k, &k, &L->c);
        else if (L->G.n) { mpq_set_si(m.val, kp_deg(&L->G), 1); rat_mul(&m, &m, &L->c); rat_sub(&k, &k, &m); }
    }
    bool zero = rat_is_zero(&k);
    rat_clear(&k); rat_clear(&m);
    return zero;
}

// the answer with its logs of polynomials in tan(H) written in sin and cos, when that removes
// every ln|cos(H)|: 2·ln|cos(x/2)| + 2·ln|tan(x/2) + 1| = ln|sin(x) + 1|; otherwise e (taken)
static Expr* tan_logs_in_sin_cos(Expr* e, const Tower* t) {
    LogTerms T = { NULL, 0 };
    collect_terms(e, &T);
    bool any = false;
    for (int i = 0; i < T.n; i++) {
        LogTerm* L = &T.t[i];
        if (!L->arg || is_cos_h(L->arg, t) || !tan_poly(L->arg, t->v, &L->G)) continue;
        if (kp_deg(&L->G) > 0) any = true;
        else lpoly_free(&L->G);                           // ln of a constant: kept as it is
    }
    Expr* out = e;
    if (any && logs_cancel(&T, t)) {
        out = create_const(0);
        bool ok = true;
        for (int i = 0; ok && i < T.n; i++) {
            const LogTerm* L = &T.t[i];
            if (L->arg && is_cos_h(L->arg, t)) continue;
            Expr* term = NULL;
            if (L->arg && L->G.n) ok = rewritten_log(L, t, &term);
            else term = expr_copy(L->e);                  // kept as it was
            if (term) out = plus(out, term);
        }
        if (ok) free_expr(e);
        else { free_expr(out); out = e; }                 // never partial: the tan form stays
    }
    log_terms_free(&T);
    return out;
}

// steps 3-4 once the normal part is reduced (g) and its logs found (L, taken)
static RischResult finish_tan(const ExpForm* F, const Reduced* R, LogPart* L) {
    KFrac rem = remainder_frac(F, R, L);
    TanRest rest;
    bool ok = tan_rest(&rem, &rest);
    kfrac_free(&rem);
    if (!ok) return failure(RISCH_UNSUPPORTED, "internal check failed: the remainder in tan has another denominator");
    Expr* acc = create_const(0);
    RischResult r = integrate_tan_rest(&rest, &F->t, &acc);
    lpoly_free(&rest.P); lpoly_free(&rest.S);
    if (r.status != RISCH_OK) { free_expr(acc); return r; }
    if (R->g.n.n) acc = plus(acc, kfrac_expr(&R->g, &F->t));
    if (L->e) { acc = plus(acc, L->e); L->e = NULL; }
    r.result = tan_logs_in_sin_cos(acc, &F->t);
    return r;
}

static RischResult integrate_tan(const ExpForm* F) {
    SpecialSplit z = split_special(&F->den);
    Reduced R = { { { 0, 0, NULL }, kp_one() }, { 0, 0, NULL }, kp_one() };
    LogPart L = { NULL, { { 0, 0, NULL }, kp_one() }, NULL };
    ResidueStatus st = RESIDUE_OK;
    if (kp_deg(&z.E) > 0) {                               // a normal part: Hermite, then residues
        LPoly s = normal_numerator(F, &z);
        reduced_free(&R);
        R = hermite_reduce(&s, &(ExpForm){ F->t, { 0, 0, NULL }, z.E, NULL });
        lpoly_free(&s);
        log_part_free(&L);
        st = residue_log_part(&(KFrac){ R.h, R.d }, &F->t, &L);
    }
    RischResult r;
    switch (st) {
        case RESIDUE_NONCONSTANT: r = (RischResult){ RISCH_NONELEMENTARY, kfrac_expr(&(KFrac){ R.h, R.d }, &F->t), RESIDUES_NOT_CONSTANT }; break;
        case RESIDUE_UNSUPPORTED: r = failure(RISCH_UNSUPPORTED, L.reason); break;
        default: r = finish_tan(F, &R, &L); break;
    }
    log_part_free(&L);
    reduced_free(&R);
    lpoly_free(&z.S); lpoly_free(&z.E);
    return r;
}

// sin/cos: exponential sums first (expsum.c); a sum in a denominator goes to θ = tan(H)
static RischResult integrate_trig(const Expr* f) {
    RischResult r = integrate_expsum(f);
    if (r.status != RISCH_UNSUPPORTED || r.reason != EXPSUM_SUM_DIVISOR) return r;
    ExpForm F = tan_form(f);
    r = F.unsupported ? failure(RISCH_UNSUPPORTED, F.unsupported) : integrate_tan(&F);
    exp_form_free(&F);
    return r;
}

// sin/cos go to exponential sums (expsum.c); so do exp-only integrands outside one tower
// (exp(x)·exp(x²), exp(x) + exp(x²)), keeping the tower's reason when that fails too
static bool use_expsum(const Expr* f, const ExpForm* F) {
    return expr_has_trig(f) || (F->unsupported && !expr_has_ln(f));
}

RischResult integrate_exp(const Expr* f) {
    if (expr_has_trig(f)) return integrate_trig(f);
    ExpForm F = exp_form(f);
    if (use_expsum(f, &F)) {
        RischResult r = integrate_expsum(f);
        if (r.status == RISCH_UNSUPPORTED) r.reason = F.unsupported;
        exp_form_free(&F);
        return r;
    }
    if (F.unsupported && expr_has_ln(f)) {                // two towers: ln over exp or ln (risch2.c)
        RischResult r2 = integrate_tower2(f);
        if (r2.reason != TOWER2_NOT_APPLICABLE) { exp_form_free(&F); return r2; }
    }
    if (F.unsupported) { const char* why = F.unsupported; exp_form_free(&F); return failure(RISCH_UNSUPPORTED, why); }
    RischResult r = kp_is_one(&F.den) ? integrate_poly(&F, f) : integrate_fraction(&F, f);
    exp_form_free(&F);
    return r;
}

// d/dx of sin/cos expressions: exponential sums, or θ = tan(B) for a sum in a denominator
// (tan(x/2) would be exact too, but unreadable: d/dx 1/(1 + cos(x)) = 1/2·tan(x/2)^3 + ...)
static Expr* differentiate_trig(const Expr* f, const char** reason) {
    Expr* r = differentiate_expsum(f, reason);
    if (r || *reason != EXPSUM_SUM_DIVISOR) return r;
    ExpForm F = tan_form(f);
    *reason = F.unsupported ? F.unsupported : F.t.half ? EXPSUM_SUM_DIVISOR : NULL;
    if (!*reason) {
        KFrac q = { F.f, F.den }, d = kfrac_derivative(&q, &F.t);
        r = kfrac_expr(&d, &F.t);
        kfrac_free(&d);
    }
    exp_form_free(&F);
    return r;
}

Expr* differentiate_exp(const Expr* f, const char** reason) {
    if (expr_has_trig(f)) return differentiate_trig(f, reason);
    ExpForm F = exp_form(f);
    Expr* r = NULL;
    const char* ignored = NULL;                           // on failure, the tower's reason is the one reported
    if (use_expsum(f, &F) && (r = differentiate_expsum(f, &ignored))) { exp_form_free(&F); return r; }
    if (F.unsupported) *reason = F.unsupported;
    else {
        KFrac q = { F.f, F.den }, d = kfrac_derivative(&q, &F.t);
        r = kfrac_expr(&d, &F.t);
        kfrac_free(&d);
    }
    exp_form_free(&F);
    return r;
}

Expr* simplify_exp(const Expr* f, const char** reason) {
    if (expr_has_trig(f)) return simplify_expsum(f, reason);
    ExpForm F = exp_form(f);
    Expr* r = NULL;
    const char* ignored = NULL;
    if (use_expsum(f, &F) && (r = simplify_expsum(f, &ignored))) { exp_form_free(&F); return r; }
    if (F.unsupported) *reason = F.unsupported;
    else r = kfrac_expr(&(KFrac){ F.f, F.den }, &F.t);
    exp_form_free(&F);
    return r;
}
