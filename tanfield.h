#ifndef TANFIELD_H
#define TANFIELD_H

#include "kpoly.h"

// Q(x)(θ), θ = tan(H) for a nonconstant polynomial H, with θ' = H'·(1 + θ²): every rational
// function of x, sin(k·B) and cos(k·B) (k integers, B one polynomial), since
//   cos(jH) + i·sin(jH) = (1 + iθ)^j·cos(H)^j  and  cos(H)² = 1/(1 + θ²).
// H = B when every sum in f has terms of one parity under B -> B + π and f is even (Bioche's
// rule: tan(x), 1/cos(x)^2 stay in tan(x)), else H = B/2 (1/(2 + cos(x)) in tan(x/2)). The result is an ExpForm with t.kind = TOWER_TAN,
// t.v = H and f/den reduced, or a reason why e is outside this setting.
ExpForm tan_form(const Expr* e);

LPoly tan_special(void);                              // 1 + θ², the special polynomial: D(1 + θ²) = 2H'·θ·(1 + θ²)

// e as a polynomial in θ = tan(H) (H given), if it is one
bool tan_poly(const Expr* e, const Polynomial* H, LPoly* p);

#endif
