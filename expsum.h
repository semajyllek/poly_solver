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
extern const char EXPSUM_SUM_DIVISOR[];               // the reason for a divisor with several terms (integrate_exp then tries tan)
Expr* differentiate_expsum(const Expr* f, const char** reason);
Expr* simplify_expsum(const Expr* f, const char** reason);
// f = c_0 + Σ_{w≠0} c_w·e^w seen from outside: c_0 (real), the product of every coefficient's
// denominator (f is finite away from its real roots; at one it may still be, as sin(x)/x),
// and per end (0: x → -∞, 1: x → +∞) whether every w ≠ 0 has Re w → -∞ there
typedef struct {
    Polynomial *c0_num, *c0_den;
    Polynomial* poles;
    bool decays[2];
    bool pure;    // f = c_0
} ExpSumShape;

bool expsum_shape(const Expr* f, ExpSumShape* out, const char** why);
void expsum_shape_free(ExpSumShape* s);
bool small_integer(const Expr* e, long* k);           // e is an integer constant with |k| <= 256

#endif
