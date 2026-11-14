#ifndef CONSTFOLD_H
#define CONSTFOLD_H

#include "expression.h"

// e with x replaced by the rational x0, folded exactly: products collect their rational
// factor (sqrt(m)·sqrt(m) = m), sums collect their rational part and merge like terms,
// ln|c·m^e| splits into logs of rationals, and atan takes its exact values at 0, ±1,
// ±sqrt(3), ±1/sqrt(3), ±(sqrt(2) ∓ 1) and ±(2 ∓ sqrt(3)). e must not contain a RootSum.
Expr* expr_at(const Expr* e, const Rational* x0);

// the folded sum of n owned constant terms. ln and atan terms over one Q(sqrt(m)) with the
// same non-rational coefficient merge (qsurd.c), every log argument is made canonical
// (> 1, content 1, no factor sqrt(m), both parts positive, not a perfect power), and
// terms print as logs, atans, the rest, then multiples of π, then the rational part.
Expr* const_fold_sum(Expr** terms, int n);

// the leading coefficient c and degree of e as x → ∞ (e ~ c·x^deg), c a folded constant; NULL
// when e is not built from x, constants and + − · / integer powers, or its leading terms cancel
Expr* expr_lead(const Expr* e, int* deg);

#endif
