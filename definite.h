#ifndef DEFINITE_H
#define DEFINITE_H

#include "expression.h"

// ∫_a^b f dx for a rational function f; a and b may be ±∞
typedef struct {
    Expr* exact;       // exact value, or NULL when only the numeric value is available
    double approx;     // adaptive quadrature, independent of the antiderivative
    const char* error; // set when there is no value at all (a pole in [a, b], not rational)
    const char* note;  // set when exact is NULL but approx is valid: the reason
} DefiniteResult;

// a bound: the rational x, or ±∞ when inf = ±1 (x is then unused)
typedef struct {
    const Rational* x;
    int inf;
} Endpoint;

typedef struct {
    Endpoint a, b;
} Bounds;

DefiniteResult integrate_definite(Expr* f, const Bounds* ab);
void definite_free(DefiniteResult* r);

#endif
