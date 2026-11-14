// The Risch algorithm for one exponential: θ = exp(v), v ∈ Q[x] nonconstant, integrand
// f = Σ a_k·θ^k with a_k ∈ Q(x). By Liouville's theorem (and θ transcendental over Q(x)),
// ∫ f is elementary iff ∫ f = ∫ a_0 + Σ_{k≠0} y_k·θ^k with y_k ∈ Q(x), where each y_k solves
// the Risch differential equation y' + k·v'·y = a_k. Here k·v' is a polynomial, so:
//   * y has poles only at poles of a_k, and ν_p(y) = ν_p(a_k) + 1 there (y' dominates).
//     A simple pole of a_k therefore has no solution: ∫ exp(x)/x is not elementary.
//   * otherwise y = q/h with h = gcd(D, D') (D = den a_k) and q a polynomial solving
//     h·q' + (k·v'·h - h')·q = N·h²/D, whose leading term fixes deg q = deg C - deg B;
//     q follows term by term from the top, and a nonzero remainder means no solution.
#include "risch.h"
#include "expfield.h"
#include "integrate.h"
#include "modgcd.h"
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

static RischResult failure(RischStatus s, const char* why) { return (RischResult){ s, NULL, why }; }

// acc + e with e's own sum spliced in term by term (so each term's sign prints as a - b); takes both
static Expr* append_sum(Expr* acc, Expr* e) {
    if (e->type == EXPR_BINOP && (e->binop.op == OP_ADD || e->binop.op == OP_SUB)) {
        Expr *l = e->binop.left, *r = e->binop.right;
        BinOp op = e->binop.op;
        e->binop.left = e->binop.right = NULL;
        free_expr(e);
        acc = append_sum(acc, l);
        return op == OP_ADD ? append_sum(acc, r) : create_sub(acc, r);
    }
    return create_add(acc, e);
}

// the constant 0, as a constant or as the zero polynomial
static bool is_zero_expr(const Expr* e) {
    return (e->type == EXPR_CONST && rat_is_zero(&e->constant)) || (e->type == EXPR_POLY && poly_is_zero(e->poly));
}

// Σ_{k≠0} y_k·θ^k with y_k' + k·v'·y_k = a_k; false (with *bad the failing index) when some
// equation has no solution
static bool exp_part(const ExpForm* F, LPoly* Y, int* bad) {
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

// coefficient of θ^k in f (zero outside its range)
static QFn coeff_at(const LPoly* f, int k) {
    int i = k - f->lo;
    return i >= 0 && i < f->n ? qfn_copy(&f->a[i]) : qfn_zero();
}

// c·q for a rational constant c
static QFn qfn_scaled(const QFn* q, const Rational* c) {
    return qfn_from(poly_scale(q->num, c), poly_copy(q->den));
}

// b_1 .. b_{n+1} from the top (b[0] unused), and g0 = a_0 - b_1·L; false at a level that fails
typedef struct {
    QFn* b;   // b[k], k = 0 .. n+1
    QFn g0;
    int n;
} LogSolution;

static bool solve_log_levels(const ExpForm* F, LogSolution* S) {
    QFn L = tower_L(&F->t);
    int n = S->n;
    QFn r_next = qfn_zero();                              // r_{k+1}
    Rational d_next;                                      // d_{k+1}, set at level k
    rat_init(&d_next);
    bool ok = true;
    for (int k = n; ok && k >= 1; k--) {
        QFn a = coeff_at(&F->f, k), rL = qfn_mul(&r_next, &L);
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
        QFn a0 = coeff_at(&F->f, 0), rL = qfn_mul(&r_next, &L);
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
static LPoly log_part_poly(const LogSolution* S) {
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
            r.result = is_zero_expr(I0) ? (free_expr(I0), eb) : is_zero_expr(eb) ? (free_expr(eb), I0) : append_sum(eb, I0);
        }
    }
    lpoly_free(&B);
    for (int k = 0; k <= S.n + 1; k++) qfn_free(&S.b[k]);
    free(S.b);
    qfn_free(&S.g0);
    return r;
}

RischResult integrate_exp(const Expr* f) {
    ExpForm F = exp_form(f);
    if (F.unsupported) { const char* why = F.unsupported; exp_form_free(&F); return failure(RISCH_UNSUPPORTED, why); }
    if (F.t.kind == TOWER_LOG) { RischResult r = integrate_log(&F, f); exp_form_free(&F); return r; }
    LPoly Y;
    int bad = 0;
    if (!exp_part(&F, &Y, &bad)) {
        RischResult r = { RISCH_NONELEMENTARY, lpoly_term_expr(&F.f.a[bad], &F.t, F.f.lo + bad),
                          "has no elementary antiderivative (Risch: y' + k·v'·y = a has no rational solution)" };
        exp_form_free(&F);
        return r;
    }
    RischResult r = { RISCH_OK, NULL, NULL };
    if (!check_exp_part(&F, &Y)) r = failure(RISCH_UNSUPPORTED, "internal check failed: the exp part does not differentiate back");
    Expr* a0 = a0_expr(&F.f);
    Expr* I0 = r.status == RISCH_OK ? integrate(a0, "x") : NULL;   // rational part (self-checked)
    if (r.status == RISCH_OK && !I0) r = failure(RISCH_UNSUPPORTED, integrate_error());
    if (r.status == RISCH_OK) {
        Expr* ey = lpoly_expr(&Y, &F.t);
        r.result = is_zero_expr(I0) ? (free_expr(I0), ey) : is_zero_expr(ey) ? (free_expr(ey), I0) : append_sum(ey, I0);
    }
    else free_expr(I0);
    free_expr(a0);
    lpoly_free(&Y);
    exp_form_free(&F);
    return r;
}

Expr* differentiate_exp(const Expr* f, const char** reason) {
    ExpForm F = exp_form(f);
    Expr* r = NULL;
    if (F.unsupported) *reason = F.unsupported;
    else {
        LPoly d = lpoly_derivative(&F.f, &F.t);
        r = lpoly_expr(&d, &F.t);
        lpoly_free(&d);
    }
    exp_form_free(&F);
    return r;
}

Expr* simplify_exp(const Expr* f, const char** reason) {
    ExpForm F = exp_form(f);
    Expr* r = NULL;
    if (F.unsupported) *reason = F.unsupported;
    else r = lpoly_expr(&F.f, &F.t);
    exp_form_free(&F);
    return r;
}
