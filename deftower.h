#ifndef DEFTOWER_H
#define DEFTOWER_H

#include "definite.h"
#include "expfield.h"

// The singular points on [lo, hi] of f = N/D in one tower (expfield.h), θ = exp(v) or
// θ = tan(H), for definite integrals (definite.c). D must have constant coefficients, and v
// or H must be linear wherever D has roots that matter. Then:
//   * θ = tan(H) is infinite at H = π/2 + kπ: f has a pole there when deg N > deg D, else F
//     may jump there (its atan(P(θ)) terms), and those points are returned;
//   * θ = r for a real root r of D is a pole of f: x = (atan(r) + kπ - h0)/h1, or for exp
//     x = (ln(r) - v0)/v1. Such x is transcendental unless r = 0 (tan, k = 0) or r = 1 (exp),
//     and N(x, r) ≠ 0 at a transcendental x (N and D are coprime, and N(·, r) is a nonzero
//     rational function with algebraic coefficients); at the rational point N is tested exactly.

void pi_50(Rational* p);                               // |π - p| < 1e-50

// lo <= hi, both finite
typedef struct {
    const Endpoint *lo, *hi;
} TowerEnds;

typedef struct {
    double* cross;        // tan: the points θ = ∞ strictly inside, ascending (f finite there)
    int n;
    bool at_lo, at_hi;    // tan: an end is such a point
} TowerPoints;

// NULL when f has no pole on [lo, hi] (out filled), else why: a divergence, or a case that
// is not supported or cannot be decided
const char* tower_points(const ExpForm* F, const TowerEnds* ends, TowerPoints* out);
void tower_points_free(TowerPoints* p);

#endif
