#ifndef RISCH_H
#define RISCH_H

#include "expression.h"

// Integration over Q(x)(exp(v)) by the Risch algorithm, for integrands Σ a_k·exp(k·v) with
// a_k in Q(x) and v a nonconstant polynomial (expfield.h). The answer is exact, or a proof
// that no elementary antiderivative exists, or "unsupported" for input outside the setting.
typedef enum { RISCH_OK, RISCH_NONELEMENTARY, RISCH_UNSUPPORTED } RischStatus;

typedef struct {
    RischStatus status;
    Expr* result;       // RISCH_OK: the antiderivative; RISCH_NONELEMENTARY: the term a·exp(k·v)
                        // whose integral is not elementary
    const char* reason; // otherwise: why (a static string)
} RischResult;

RischResult integrate_exp(const Expr* f);

// d/dx of f in the same setting; NULL (and *reason) when f is outside it
Expr* differentiate_exp(const Expr* f, const char** reason);

// f in canonical form Σ a_k·exp(k·v); NULL (and *reason) when f is outside the setting
Expr* simplify_exp(const Expr* f, const char** reason);

#endif
