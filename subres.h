#ifndef SUBRES_H
#define SUBRES_H

#include "polynomial.h"

// Bivariate polynomial Σ c[k](t) · x^k, each c[k] a univariate polynomial in t.
typedef struct {
    Polynomial** c;
    int deg; // degree in x (0 for the zero polynomial)
} PolyTX;

PolyTX ptx_alloc(int deg); // all-zero coefficients
PolyTX ptx_copy(const PolyTX* a);
void   ptx_free(PolyTX* a);
bool   ptx_is_zero(const PolyTX* a);

// The Rothstein–Trager pair for ∫ a/d: d(x) and g(t, x) = a(x) - t·d'(x).
// Requires deg a < deg d and deg d >= 2.
typedef struct {
    const Polynomial* d;
    const Polynomial* a;
} RTPair;

PolyTX rt_g(const RTPair* p); // a - t·d' as a bivariate polynomial

// The full subresultant chain of d and g in x over Q[t]: S[j] for j = 0..m (m = deg d),
// S[m] = d, S[m-1] = g, S[j] the Collins subresultant of index j (zero inside gaps).
// The resultant R(t) = res_x(d, g) is S[0].c[0].
typedef struct {
    PolyTX* S;
    int m;
} SubresChain;

SubresChain rt_chain(const RTPair* p);
void        chain_free(SubresChain* ch);

#endif
