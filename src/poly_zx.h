#ifndef POLY_ZX_H
#define POLY_ZX_H

#include <gmp.h>
#include <stdbool.h>
#include "polynomial.h"
#include "poly_zp.h"

// Dense univariate polynomials over the integers with GMP coefficients.
// c[i] is the coefficient of x^i; deg == -1 is the zero polynomial.
typedef struct {
    mpz_t* c;
    int deg;
} ZxPoly;

typedef struct {
    ZxPoly* items;
    int count;
} ZxFactorList;

// lifecycle and conversion
ZxPoly zx_alloc(int deg); // all-zero coefficients up to deg
ZxPoly zx_copy(const ZxPoly* a);
void   zx_free(ZxPoly* a);
void   zx_trim(ZxPoly* a);
ZxPoly zx_from_poly(const Polynomial* p); // p must have integer coefficients
Polynomial* zx_to_poly(const ZxPoly* a);

// arithmetic
ZxPoly zx_mul(const ZxPoly* a, const ZxPoly* b);
bool   zx_divexact(const ZxPoly* a, const ZxPoly* b, ZxPoly* q); // true iff b | a in Z[x]
void   zx_make_primitive(ZxPoly* a);              // divide by content, make lc > 0
void   zx_mods(ZxPoly* a, const mpz_t m);         // reduce into the symmetric range (-m/2, m/2]
ZpPoly zx_to_zp(const ZxPoly* a, const ZpField* F);
ZxPoly zx_from_zp(const ZpPoly* a);               // coefficients in [0, p)

// Berlekamp–Zassenhaus: irreducible factors over Z of a primitive,
// square-free f with positive leading coefficient
ZxFactorList zx_factor_squarefree(const ZxPoly* f);
void zx_factor_list_free(ZxFactorList* L);

#endif
