#ifndef SUBRES_H
#define SUBRES_H

#include "polynomial.h"

// Bivariate polynomial Σ c[k](t) · x^k, each c[k] a univariate polynomial in t.
typedef struct {
    Polynomial** c;
    int deg; // degree in x
} PolyTX;

PolyTX ptx_alloc(int deg); // all-zero coefficients
PolyTX ptx_copy(const PolyTX* a);
void   ptx_free(PolyTX* a);

// The Rothstein–Trager pair for ∫ a/d: d(x) and g(t, x) = a(x) - t·d'(x).
// Requires deg a < deg d and deg d >= 2.
typedef struct {
    const Polynomial* d;
    const Polynomial* a;
} RTPair;

PolyTX      rt_g(const RTPair* p);                 // a - t·d' as a bivariate polynomial
Polynomial* rt_resultant(const RTPair* p);          // R(t) = res_x(d, a - t·d')
PolyTX      rt_subresultant(const RTPair* p, int j); // S_j(t, x), 1 <= j < deg d - 1

// polynomial of degree < n through (ts[i], vs[i]), exact Newton interpolation
Polynomial* interpolate(const Rational* ts, const Rational* vs, int n);

#endif
