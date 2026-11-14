#ifndef QSURD_H
#define QSURD_H

#include "expression.h"

// p + q·sqrt(m) in Q(sqrt(m)): exact arithmetic for folding constants such as
// ln|2 + sqrt(2)| - ln|2 - sqrt(2)| = ln(3 + 2·sqrt(2))
typedef struct {
    Rational p, q;
    long m; // the radicand, 0 when q = 0
} QSurd;

void qs_init(QSurd* a);
void qs_clear(QSurd* a);
void qs_set(QSurd* r, const QSurd* a);

// e as p + q·sqrt(m): a rational, c·sqrt(m) or c·sqrt(m)/d, and sums of these with one m
bool qs_parse(const Expr* e, QSurd* out);

// a ∘ b; false (r unchanged) when a and b lie in different fields, or b = 0 for division
bool qs_add(QSurd* r, const QSurd* a, const QSurd* b);
bool qs_sub(QSurd* r, const QSurd* a, const QSurd* b);
bool qs_mul(QSurd* r, const QSurd* a, const QSurd* b);
bool qs_div(QSurd* r, const QSurd* a, const QSurd* b);

int qs_sign(const QSurd* a);           // exact
bool qs_sqrt(QSurd* h, const QSurd* g); // h > 0 with h² = g, both with integer parts; false if none
Expr* qs_expr(const QSurd* a);         // p + q·sqrt(m), rational part first

#endif
