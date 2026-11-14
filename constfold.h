#ifndef CONSTFOLD_H
#define CONSTFOLD_H

#include "expression.h"

// e with x replaced by the rational x0, folded exactly: products collect their rational
// factor (sqrt(m)·sqrt(m) = m), sums collect their rational part and merge like terms,
// ln|c·m^e| splits into logs of rationals, and atan takes its exact values at 0, ±1,
// ±sqrt(3) and ±1/sqrt(3). e must not contain a RootSum.
Expr* expr_at(const Expr* e, const Rational* x0);

// the folded sum of n owned constant terms: other terms first, then c·ln(v) (like v merged),
// then the multiple of π, then the rational part
Expr* const_fold_sum(Expr** terms, int n);

// the leading coefficient c and degree of e as x → ∞ (e ~ c·x^deg), c a folded constant; NULL
// when e is not built from x, constants and + − · / integer powers, or its leading terms cancel
Expr* expr_lead(const Expr* e, int* deg);

#endif
