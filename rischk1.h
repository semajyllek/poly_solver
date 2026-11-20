#ifndef RISCHK1_H
#define RISCHK1_H

#include "risch.h"
#include "kpoly.h"

// The one-tower Risch machinery of risch.c over K1 = Q(x)(θ1), for the two-level integrator
// (risch2.c): f is an ExpForm (expfield.h) in the tower F->t.

RischResult failure(RischStatus s, const char* why);
bool is_zero_expr(const Expr* e);
Expr* plus(Expr* acc, Expr* e);                       // acc + e, either possibly 0; takes both

// θ = exp(v): y_k with y_k' + k·v'·y_k = a_k for every k ≠ 0; false (*bad the index) when one fails
bool exp_part(const ExpForm* F, LPoly* Y, int* bad);

// θ = ln(u), f a polynomial in θ: b_1 .. b_{n+1} from the top (b[0] unused) and g0 = a_0 - b_1·L;
// false at a level that needs a logarithm other than θ
typedef struct {
    QFn* b;   // b[k], k = 0 .. n+1
    QFn g0;
    int n;
} LogSolution;
bool solve_log_levels(const ExpForm* F, LogSolution* S);
LPoly log_part_poly(const LogSolution* S);            // Σ_{k>=1} b_k θ^k

// Hermite reduction of b/F->den: ∫ b/den = g + ∫ h/d, d square-free, deg h < deg d
typedef struct {
    KFrac g;     // the rational part of the integral
    LPoly h, d;  // the remaining h/d
} Reduced;
LPoly proper_part(const ExpForm* F);                  // b with f = p + b/den, deg b < deg den
Reduced hermite_reduce(const LPoly* b, const ExpForm* F);
void reduced_free(Reduced* R);

// ∫ f over K1 for a polynomial (Laurent for exp) f, or a fraction; f the integrand's expression
RischResult integrate_poly(const ExpForm* F, const Expr* f);
RischResult integrate_fraction(const ExpForm* F, const Expr* f);

#endif
