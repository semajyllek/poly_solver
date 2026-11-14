#ifndef SURD_H
#define SURD_H

#include "expression.h"

// Exact real constants built from rationals, radicals and trig values at rational
// multiples of π: one monomial is c · m^e · sqrt(s) · trig(angle·π), where
//   m >= 1 is not a perfect power and 0 <= e < 1 (m = 1, e = 0 for none),
//   sqrt(s), s squarefree, is used only beside a radical that is not a square root
//   (2^(1/5)·sqrt(5)); otherwise square roots merge into m^e (sqrt(2)·sqrt(3) = sqrt(6)),
//   trig is symbolic only off the table, with angle in (0, 1/2).
typedef struct {
    Rational c;
    mpz_t m;
    Rational e;
    int s;
    TrigKind kind;
    Rational angle;
} Mono;

// R^e for R > 0 and any rational e
typedef struct {
    const Rational* R;
    const Rational* e;
} Power;

// x = c · R^e, with the radical pulled into canonical form
void mono_init_radical(Mono* x, const Rational* c, const Power* pw);
// x <- x · trig(angle·π); at most one trig factor per monomial
void mono_mul_trig(Mono* x, TrigKind kind, const Rational* angle);
void mono_clear(Mono* x);

// Σ monomials with like terms merged and zeros dropped
typedef struct {
    Mono* t;
    int count;
} RealConst;

RealConst rc_zero(void);
void rc_add(RealConst* k, Mono* x); // takes ownership of x's contents
bool rc_is_zero(const RealConst* k);
void rc_clear(RealConst* k);

// k · body (a parenthesized sum when k has several monomials); takes body; NULL when k = 0
Expr* rc_scale(const RealConst* k, Expr* body);
Expr* rc_x_pow(int power); // the polynomial x^power
// acc + Σ_i mono_i · x^power, one term per monomial (power < 0: the bare constant)
Expr* rc_add_terms(Expr* acc, const RealConst* k, int power);

#endif
