#ifndef ROOTS_H
#define ROOTS_H

#include "polynomial.h"

// An isolating interval: exactly one real root lies in the open interval (lo, hi),
// or the root is exactly lo when lo == hi. Endpoints are exact rationals.
typedef struct {
    Rational lo, hi;
} RootInterval;

typedef struct {
    RootInterval* items; // sorted ascending
    int count;
} RootIntervalList;

// Isolate all real roots of a square-free polynomial with the Descartes rule of
// signs and bisection (Vincent–Collins–Akritas). Exact arithmetic throughout.
RootIntervalList isolate_real_roots(const Polynomial* f);

// Shrink an isolating interval by exact bisection until hi - lo <= width.
void refine_root(RootInterval* r, const Polynomial* f, const Rational* width);

void root_interval_list_free(RootIntervalList* L);

#endif
