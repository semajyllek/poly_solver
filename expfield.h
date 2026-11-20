#ifndef EXPFIELD_H
#define EXPFIELD_H

#include "expression.h"

// The differential field Q(x)(θ) for one transcendental θ over Q(x), either
//   θ = exp(v), v a nonconstant polynomial (θ' = v'·θ), with Laurent polynomials Σ a_k·θ^k, or
//   θ = ln(u), u a nonconstant rational function (θ' = u'/u), with polynomials Σ_{k>=0} a_k·θ^k,
// coefficients a_k in Q(x): x²·exp(3x), exp(x)/x, x·ln(x)^2, ln(x^2 + 1)/x, ...; and quotients of
// these (1/(1 + exp(x)), 1/(x·ln(x))), kept as reduced fractions over K[θ] (kpoly.h)

// a rational function num/den in lowest terms with a monic den
typedef struct {
    Polynomial *num, *den;
} QFn;

// Σ a[i]·θ^(lo + i)
typedef struct {
    int lo, n;
    QFn* a;
} LPoly;

typedef enum { TOWER_EXP, TOWER_LOG, TOWER_TAN } TowerKind;

// θ = exp(v), θ = ln(u) with L = u'/u, or θ = tan(v) (θ' = v'·(1 + θ²), polynomials in θ; tanfield.h)
typedef struct {
    TowerKind kind;
    Polynomial* v;   // TOWER_EXP, TOWER_TAN
    Polynomial *un, *ud, *Ln, *Ld; // TOWER_LOG: u = un/ud and L = Ln/Ld in lowest terms
    bool half;       // TOWER_TAN: v = B/2 for the base B of the arguments (tanfield.h)
} Tower;

// the tower of e and e as a Laurent polynomial in it, or a reason (static string) why e is
// outside this setting (division by a sum with exp, exp of a non-polynomial, ln, ...)
typedef struct {
    Tower t;      // t.v is the constant 1 when e has no exp of a nonconstant argument
    LPoly f, den; // e = f/den in lowest terms, den monic in θ (θ ∤ den for exp); den = 1 mostly
    const char* unsupported;
} ExpForm;

ExpForm exp_form(const Expr* e);
void exp_form_free(ExpForm* F);
QFn tower_L(const Tower* t);                          // u'/u for TOWER_LOG (a copy)

QFn qfn_zero(void);
QFn qfn_copy(const QFn* a);
void qfn_free(QFn* a);
bool qfn_is_zero(const QFn* a);
Expr* qfn_expr(const QFn* a);                         // a as EXPR_POLY or EXPR_RATIONAL_FN
QFn qfn_from(Polynomial* num, Polynomial* den);       // takes both, reduces
QFn qfn_poly(const Polynomial* p);                    // p/1 (copies)
QFn qfn_add(const QFn* a, const QFn* b, int sign);    // a ± b
QFn qfn_mul(const QFn* a, const QFn* b);
LPoly lpoly_term(QFn a, int k);                       // a·θ^k; takes a
LPoly lpoly_add(const LPoly* f, const LPoly* g, int sign);
LPoly lpoly_mul(const LPoly* f, const LPoly* g);
QFn qfn_inv(const QFn* a);                            // 1/a for a ≠ 0
QFn qfn_derivative(const QFn* a);                     // d/dx

LPoly lpoly_derivative(const LPoly* f, const Tower* t); // (a·θ^k)' = (a' + k·v'·a)·θ^k
bool lpoly_equal(const LPoly* f, const LPoly* g);
void lpoly_free(LPoly* f);
Expr* lpoly_expr(const LPoly* f, const Tower* t);     // Σ a_k·exp(k·v), 0 when empty
Expr* lpoly_expr_order(const LPoly* f, const Tower* t, bool down);
Expr* lpoly_term_expr(const QFn* a, const Tower* t, int k); // one term a·exp(k·v)

#endif
