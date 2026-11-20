#ifndef DEFINITE_H
#define DEFINITE_H

#include "expression.h"

// ∫_a^b f dx for a rational function f, or an exponential sum (exp, sin, cos with rational
// coefficients, no sum in a denominator); a and b may be ±∞ or rational multiples of π
typedef struct {
    Expr* exact;       // exact value, or NULL when only the numeric value is available
    double approx;     // adaptive quadrature, independent of the antiderivative
    const char* error; // set when there is no value at all (a pole in [a, b], not rational)
    const char* note;  // set when exact is NULL but approx is valid: the reason (or, with exact, why
                       // approx could not be compared: a value beyond double precision)
} DefiniteResult;

// a bound: the rational x (x·π when pi, x ≠ 0), or ±∞ when inf = ±1 (x is then unused)
typedef struct {
    const Rational* x;
    int inf;
    bool pi;
} Endpoint;

typedef struct {
    Endpoint a, b;
} Bounds;

DefiniteResult integrate_definite(Expr* f, const Bounds* ab);

// a bound from text: inf, +inf, oo, -inf, -oo; q·π as pi, -pi/2, 3pi/4, 2*pi; or a constant
// expression ("1/2", "-3") evaluated exactly. The rational goes to *out, which end->x should point to
bool parse_endpoint(Rational* out, Endpoint* end, const char* text);
// an integrand from text: the raw tree when it has exp, sin, cos or ln, else simplified
Expr* parse_integrand(const char* text);
void definite_free(DefiniteResult* r);

#endif
