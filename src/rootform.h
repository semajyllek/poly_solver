#ifndef ROOTFORM_H
#define ROOTFORM_H

#include "expression.h"

// The roots R^(1/n)·e^(iπ·j/L), j in js, of an irreducible P whose roots are known exactly:
// a binomial p·x^n + q, or Q(x^s) with Q = q·b^m·Φ_N(y/b) for rational b > 0.
typedef struct {
    Rational R;
    int n;
    int L;
    int* js;
    int count;
} RootSet;

bool  known_roots(const Polynomial* P, RootSet* rs); // false (rs untouched) for any other P
void  root_set_free(RootSet* rs);
// ρ·cos(πθ) (TRIG_COS) or ρ·sin(πθ) (TRIG_SIN) for ρ = R^(1/n), exact; NULL when zero
Expr* root_part(const RootSet* rs, const Rational* theta, TrigKind kind);

// ∫ a/P in real radical/trig form when the roots of the irreducible P (deg P >= 3,
// deg a < deg P) are known exactly: P a binomial p·x^n + q, or a rational rescaling
// b^n·Φ_N(x/b) of a cyclotomic polynomial. NULL otherwise.
Expr* root_log_part(const Polynomial* a, const Polynomial* P);

#endif
