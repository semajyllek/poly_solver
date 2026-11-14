#ifndef ROOTFORM_H
#define ROOTFORM_H

#include "expression.h"

// ∫ a/P in real radical/trig form when the roots of the irreducible P (deg P >= 3,
// deg a < deg P) are known exactly: P a binomial p·x^n + q, or a rational rescaling
// b^n·Φ_N(x/b) of a cyclotomic polynomial. NULL otherwise.
Expr* root_log_part(const Polynomial* a, const Polynomial* P);

#endif
