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

RischResult integrate_exp(const Expr* f) {
    ExpForm F = exp_form(f);
    if (F.unsupported) { const char* why = F.unsupported; exp_form_free(&F); return failure(RISCH_UNSUPPORTED, why); }
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
        r.result = is_zero_expr(I0) ? (free_expr(I0), ey) : is_zero_expr(ey) ? (free_expr(ey), I0) : create_add(ey, I0);
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
