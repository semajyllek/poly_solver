#ifndef RISCH2_H
#define RISCH2_H

#include "risch.h"

// Two levels: θ2 = ln(u) over K1 = Q(x)(θ1), θ1 = exp(v) or ln(u1), and u ∈ K1, for integrands
// that are polynomials in θ2 with coefficients in K1: exp(x)·ln(x), ln(ln(x))/x, ln(x)·ln(x+1).
// The answer is exact (checked), or a proof that none is elementary, or "unsupported".
// The reason is TOWER2_NOT_APPLICABLE when f does not have this shape at all.
extern const char TOWER2_NOT_APPLICABLE[];

RischResult integrate_tower2(const Expr* f);

#endif
