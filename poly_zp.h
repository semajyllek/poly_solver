#ifndef POLY_ZP_H
#define POLY_ZP_H

#include <gmp.h>
#include <stdbool.h>

// Dense univariate polynomials over the prime field F_p, p < 2^31.
// c[i] is the coefficient of x^i; deg == -1 is the zero polynomial.
typedef struct {
    unsigned long* c;
    int deg;
} ZpPoly;

typedef struct {
    unsigned long p;
} ZpField;

// arithmetic modulo a fixed polynomial (used by powmod)
typedef struct {
    const ZpField* F;
    const ZpPoly* mod;
} ZpRing;

typedef struct {
    ZpPoly s, t; // s*a + t*b == 1
} ZpBezout;

typedef struct {
    ZpPoly* items;
    int count;
} ZpFactorList;

// lifecycle
ZpPoly zp_alloc(int deg); // all-zero coefficients up to deg
ZpPoly zp_copy(const ZpPoly* a);
void   zp_free(ZpPoly* a);
void   zp_trim(ZpPoly* a);

// field and ring arithmetic
unsigned long zp_inv(unsigned long a, const ZpField* F);
ZpPoly zp_add(const ZpPoly* a, const ZpPoly* b, const ZpField* F);
ZpPoly zp_sub(const ZpPoly* a, const ZpPoly* b, const ZpField* F);
ZpPoly zp_mul(const ZpPoly* a, const ZpPoly* b, const ZpField* F);
ZpPoly zp_rem(const ZpPoly* a, const ZpPoly* b, const ZpField* F);
ZpPoly zp_quo(const ZpPoly* a, const ZpPoly* b, const ZpField* F);
ZpPoly zp_gcd(const ZpPoly* a, const ZpPoly* b, const ZpField* F); // monic
ZpPoly zp_deriv(const ZpPoly* a, const ZpField* F);
void   zp_make_monic(ZpPoly* a, const ZpField* F);
ZpPoly zp_powmod(const ZpPoly* base, const mpz_t e, const ZpRing* R);
ZpBezout zp_xgcd(const ZpPoly* a, const ZpPoly* b, const ZpField* F); // requires gcd(a,b)=1

// complete factorization of a square-free polynomial into monic irreducibles
// (distinct-degree + Cantor–Zassenhaus equal-degree splitting; p must be odd)
ZpFactorList zp_factor_squarefree(const ZpPoly* f, const ZpField* F);
void zp_factor_list_free(ZpFactorList* L);

#endif
