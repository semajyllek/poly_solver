#ifndef K2POLY_H
#define K2POLY_H

#include "kpoly.h"

// Polynomials in a second transcendental θ2 over K1 = Q(x)(θ1): coefficients are reduced
// fractions in θ1 (KFrac, kpoly.h) of one tower kind. The algorithms are kpoly_tmpl.h's, the
// same as kpoly.h's over Q(x); fractions K2Frac are reduced with a monic denominator in θ2.

typedef struct {
    KFrac* a;         // a[k]: the coefficient of θ2^k
    int n;
    TowerKind kind;   // of K1, for the arithmetic of the coefficients
} K2Poly;

KFrac kf_zero(void);
bool kf_is_zero(const KFrac* k);
KFrac kf_copy(const KFrac* k);
KFrac kf_scaled(const KFrac* k, const Rational* c, TowerKind kind);   // c·k

K2Poly k2_zeros(int n, TowerKind kind);                // n zero coefficients
K2Poly k2_const(KFrac c, TowerKind kind);              // takes c
void k2_free(K2Poly* p);
void k2_trim(K2Poly* p);
void k2_accumulate(K2Poly* p, int k, KFrac c);         // a[k] += c (taken)
int k2_deg(const K2Poly* p);                           // -1 for 0
K2Poly k2_add(const K2Poly* x, const K2Poly* y, int sign);
K2Poly k2_mul(const K2Poly* x, const K2Poly* y);
K2Poly k2_copy(const K2Poly* a);
K2Poly k2_monic(const K2Poly* a);
K2Poly k2_dtheta(const K2Poly* a);                     // d/dθ2

typedef struct {
    K2Poly q, r;
} K2Div;

typedef struct {
    K2Poly s, t;
} K2Bezout;

typedef struct {
    K2Poly* f;
    int m;
} K2SquareFree;

K2Div k2_divmod(const K2Poly* a, const K2Poly* b);     // b ≠ 0
K2Poly k2_exact_div(const K2Poly* a, const K2Poly* b);
K2Poly k2_gcd(const K2Poly* a, const K2Poly* b);       // monic
K2Bezout k2_bezout(const K2Poly* a, const K2Poly* b, const K2Poly* c);  // s·a + t·b = c, gcd(a, b) = 1
K2SquareFree k2_square_free(const K2Poly* d);          // d monic, nonconstant
void k2div_free(K2Div* d);
void k2bezout_free(K2Bezout* z);
void k2sf_free(K2SquareFree* s);

// n/d with gcd 1 and d monic in θ2
typedef struct {
    K2Poly n, d;
} K2Frac;

K2Frac k2frac_make(K2Poly n, K2Poly d);                // takes both; d ≠ 0
K2Frac k2frac_add(const K2Frac* a, const K2Frac* b, int sign);
K2Frac k2frac_mul(const K2Frac* a, const K2Frac* b);
K2Frac k2frac_inv(const K2Frac* a);                    // a ≠ 0
void k2frac_free(K2Frac* f);

#endif
