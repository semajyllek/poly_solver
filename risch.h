#ifndef RISCH_H
#define RISCH_H

#include "expression.h"

// Integration over Q(x)(θ) by the Risch algorithm, θ = exp(v) (v a nonconstant polynomial),
// θ = ln(u) (u a nonconstant rational function) or θ = tan(H) (tanfield.h), for rational
// functions of θ with coefficients in Q(x) (expfield.h); sin/cos without a sum in a
// denominator, and several exponentials, as exponential sums over Q(i)(x) (expsum.h). The answer is exact, or a proof that no elementary
// antiderivative exists, or "unsupported" for input outside the setting.
typedef enum { RISCH_OK, RISCH_NONELEMENTARY, RISCH_UNSUPPORTED } RischStatus;

typedef struct {
    RischStatus status;
    Expr* result;       // RISCH_OK: the antiderivative; RISCH_NONELEMENTARY: the part (a·exp(k·v),
                        // a reduced fraction h/d, or f) whose integral is not elementary
    const char* reason; // otherwise: why (a static string)
} RischResult;

RischResult integrate_exp(const Expr* f);

// d/dx of f in the same setting; NULL (and *reason) when f is outside it
Expr* differentiate_exp(const Expr* f, const char** reason);

// f in canonical form (Σ a_k·exp(k·v), or a reduced fraction in θ); NULL (and *reason) outside the setting
Expr* simplify_exp(const Expr* f, const char** reason);

#endif
