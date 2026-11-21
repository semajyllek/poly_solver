#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "polynomial.h"

// num/den = (rnum/rden)' + a/dstar: the rational part of the integral (Hermite reduction plus
// the integrated polynomial part) and the remaining log-part integrand, with dstar square-free
// and deg a < deg dstar. False when Hermite's exact self-check fails.
typedef struct {
    Polynomial *rnum, *rden, *a, *dstar;
} HermiteSplit;

bool hermite_split(const Polynomial* num, const Polynomial* den, HermiteSplit* out);
void hermite_split_free(HermiteSplit* h);

#endif
