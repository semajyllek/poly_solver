// One square root: expressions in x and y = sqrt(U), U a polynomial, as elements a + b·y of
// Q(x)(y) (a, b ∈ Q(x)). simplify and differentiate work for any U; integrate needs deg U <= 2:
// the conic y² = U is parametrized rationally (Euler's substitutions, or a line through a
// rational point), the integral becomes a rational one in t, and t = T(x, y) is substituted
// back. deg U >= 3 (elliptic and beyond) and conics without a rational point are unsupported.
#ifndef ALGSQRT_H
#define ALGSQRT_H

#include "expression.h"
#include "risch.h"

bool expr_has_sqrt(const Expr* e);                          // sqrt(...) or a half-integer power
RischResult integrate_sqrt(const Expr* f);                  // RISCH_OK or RISCH_UNSUPPORTED (never a proof)
Expr* differentiate_sqrt(const Expr* f, const char** why);  // NULL with *why on failure
Expr* simplify_sqrt(const Expr* f, const char** why);

#endif
