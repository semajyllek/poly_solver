// The logarithmic part of ∫ r/(d·sqrt(U)) for deg U = 2 (sqrtlog.c), d irreducible over Q,
// monic and coprime to U, deg r < deg d: Σ over the roots ρ of d of c·∫ dx/((x - ρ)·sqrt(U)).
#ifndef SQRTLOG_H
#define SQRTLOG_H

#include "expression.h"

typedef struct {
    const Polynomial *U, *r, *d;
} SqrtLog;

// the terms (real ln/atan with surds, or RootSums), or NULL with *why
Expr* sqrt_log_part(const SqrtLog* in, const char** why);

#endif
