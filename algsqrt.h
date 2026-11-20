// One square root: expressions in x and y = sqrt(U), U a polynomial, as elements a + b·y of
// Q(x)(y) (a, b ∈ Q(x)). simplify and differentiate work for any U; integrate needs deg U <= 2:
// the conic y² = U is parametrized rationally (Euler's substitutions, or a line through a
// rational point), the integral becomes a rational one in t, and t = T(x, y) is substituted
// back. A conic without a rational point takes N/(U^m·sqrt(U)) by reduction to ∫1/sqrt(U)
// (atan or ln); other forms there, and deg U >= 3 (elliptic and beyond), are unsupported.
#ifndef ALGSQRT_H
#define ALGSQRT_H

#include "expression.h"
#include "risch.h"
#include "definite.h"

bool expr_has_sqrt(const Expr* e);                          // sqrt(...) or a half-integer power
RischResult integrate_sqrt(const Expr* f);                  // RISCH_OK or RISCH_UNSUPPORTED (never a proof)
Expr* differentiate_sqrt(const Expr* f, const char** why);  // NULL with *why on failure
Expr* simplify_sqrt(const Expr* f, const char** why);

// ∫_a^b f for finite rational bounds (sqrtdef.c): exact when the antiderivative evaluates at both
// ends, always with tanh-sinh quadrature; an error for a range where sqrt(U) is not real or f has a pole
DefiniteResult integrate_definite_sqrt(const Expr* f, const Bounds* ab);

// for definite integrals (sqrtdef.c): f = a + b·sqrt(U), U NULL when f has no root, with
// a = num_a/den_a and b = num_b/den_b (denominators monic)
typedef struct {
    Polynomial *U, *num_a, *den_a, *num_b, *den_b;
} SqrtShape;
bool sqrt_shape(const Expr* f, SqrtShape* s, const char** why);
void sqrt_shape_free(SqrtShape* s);
// the value at x0 of g, an expression in x and sqrt(U), as a constant r1 + r2·sqrt(m); a limit
// through the conjugate where a denominator vanishes (((y - 1)/x at 0); NULL when g is infinite
// there or not of that form
Expr* sqrt_value_at(const Expr* g, const Polynomial* U, const Rational* x0);

#endif
