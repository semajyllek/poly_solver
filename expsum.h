#ifndef EXPSUM_H
#define EXPSUM_H

#include "risch.h"

// Exponential sums Σ c_w·e^w with c_w ∈ Q(i)(x) and exponents w = E + i·F, E, F ∈ Q[x]: every
// integrand built from +, -, ·, integer powers, division by single terms, and exp, sin, cos of
// polynomials (sin u = (e^{iu} - e^{-iu})/2i, cos u = (e^{iu} + e^{-iu})/2; products of exponentials
// add their exponents, so exp(x)·exp(x²) is one term). Since D(c·e^w) = (c' + w'·c)·e^w acts on each
// term alone, Liouville's theorem gives: ∫ is elementary iff y' + w'·y = c_w has a solution in
// Q(i)(x) for every w ≠ 0. The answer is written back with real exp, sin and cos.

RischResult integrate_expsum(const Expr* f);
Expr* differentiate_expsum(const Expr* f, const char** reason);
Expr* simplify_expsum(const Expr* f, const char** reason);

#endif
