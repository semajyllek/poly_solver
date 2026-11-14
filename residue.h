#ifndef RESIDUE_H
#define RESIDUE_H

#include "kpoly.h"

// The logarithmic part of ∫ h/d over K(θ), K = Q(x), for d monic, square-free and normal
// (gcd(d, D(d)) = 1) and deg h < deg d: the Rothstein–Trager residue test, then
//   Σ over the residues α of α·ln(gcd(d, h - α·D(d)))
// in real form: c·ln|G| for rational α; for α = a ± b·√m, a·ln|N| plus b·√m·ln|(P + √m·Q)/(P - √m·Q)|
// (real) or Rioboo's arctangents of polynomials in θ (complex, as in atan(exp(x))).
typedef enum { RESIDUE_OK, RESIDUE_NONCONSTANT, RESIDUE_UNSUPPORTED } ResidueStatus;

typedef struct {
    Expr* e;     // the terms (NULL when h = 0)
    KFrac deriv; // their derivative D(e), exactly, for the caller's remainder check
    const char* reason; // RESIDUE_UNSUPPORTED: why (a static string)
} LogPart;

// hd = h/d (not necessarily in lowest terms)
ResidueStatus residue_log_part(const KFrac* hd, const Tower* t, LogPart* out);
void log_part_free(LogPart* L);

#endif
