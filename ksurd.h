#ifndef KSURD_H
#define KSURD_H

#include "kpoly.h"

// Polynomials in θ over K(√M), K = Q(x), M a rational that is not a square (M < 0 for
// complex residues): a = p + √M·q with p, q in K[θ]. Used for g_α = gcd(d, h - α·D(d)) when
// the residue α = a + b·√M is a quadratic irrationality.
typedef struct {
    LPoly p, q;
} SPoly;

// u + √M·w
typedef struct {
    QFn u, w;
} SElt;

SPoly sp_from(const LPoly* p, const LPoly* q);        // copies both
void sp_free(SPoly* a);
int sp_deg(const SPoly* a);                            // -1 for 0
SPoly sp_gcd(const SPoly* a, const SPoly* b, const Rational* M); // monic; a or b nonzero

#endif
