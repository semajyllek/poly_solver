#ifndef EXPFIELD_H
#define EXPFIELD_H

#include "expression.h"

// The differential field Q(x)(θ) with θ = exp(v), v a nonconstant polynomial, restricted to
// Laurent polynomials Σ a_k·θ^k with coefficients a_k in Q(x): enough for integrands such as
// x²·exp(3x), exp(x)/x or (x - 1)·exp(x^2) + 1/x. (θ' = v'·θ, and θ is transcendental over Q(x).)

// a rational function num/den in lowest terms with a monic den
typedef struct {
    Polynomial *num, *den;
} QFn;

// Σ a[i]·θ^(lo + i)
typedef struct {
    int lo, n;
    QFn* a;
} LPoly;

// θ = exp(v)
typedef struct {
    Polynomial* v;
} Tower;

// the tower of e and e as a Laurent polynomial in it, or a reason (static string) why e is
// outside this setting (division by a sum with exp, exp of a non-polynomial, ln, ...)
typedef struct {
    Tower t;      // t.v is the constant 1 when e has no exp of a nonconstant argument
    LPoly f;
    const char* unsupported;
} ExpForm;

ExpForm exp_form(const Expr* e);
void exp_form_free(ExpForm* F);

QFn qfn_zero(void);
QFn qfn_copy(const QFn* a);
void qfn_free(QFn* a);
bool qfn_is_zero(const QFn* a);
Expr* qfn_expr(const QFn* a);                         // a as EXPR_POLY or EXPR_RATIONAL_FN
QFn qfn_from(Polynomial* num, Polynomial* den);       // takes both, reduces

LPoly lpoly_derivative(const LPoly* f, const Tower* t); // (a·θ^k)' = (a' + k·v'·a)·θ^k
bool lpoly_equal(const LPoly* f, const LPoly* g);
void lpoly_free(LPoly* f);
Expr* lpoly_expr(const LPoly* f, const Tower* t);     // Σ a_k·exp(k·v), 0 when empty
Expr* lpoly_term_expr(const QFn* a, const Tower* t, int k); // one term a·exp(k·v)

#endif
