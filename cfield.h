#ifndef CFIELD_H
#define CFIELD_H

#include "polynomial.h"

// Q(i)[x] and Q(i)(x): the coefficient field once sin and cos are written as complex
// exponentials (sin u = (e^{iu} - e^{-iu})/2i).

// re + i·im
typedef struct {
    Polynomial *re, *im;
} CPoly;

// (re + i·im)/den with den real and monic, gcd(re, im, den) = 1
typedef struct {
    CPoly n;
    Polynomial* den;
} CFn;

CPoly cp_zero(void);
CPoly cp_copy(const CPoly* a);
void cp_free(CPoly* a);
bool cp_is_zero(const CPoly* a);
int cp_deg(const CPoly* a);                            // -1 for 0
bool cp_equal(const CPoly* a, const CPoly* b);
CPoly cp_add(const CPoly* a, const CPoly* b, int sign);
CPoly cp_mul(const CPoly* a, const CPoly* b);
CPoly cp_neg(const CPoly* a);
CPoly cp_derivative(const CPoly* a);

CFn cfn_from(CPoly n, Polynomial* den);              // takes both, reduces; den ≠ 0
CFn cfn_real(const Polynomial* num, const Polynomial* den); // num/den (copies)
CFn cfn_const(long re_num, long im_num, long den);    // (re + i·im)/den
CFn cfn_poly(const CPoly* a);                         // a/1 (copies)
CFn cfn_zero(void);
CFn cfn_copy(const CFn* a);
void cfn_free(CFn* a);
bool cfn_is_zero(const CFn* a);
CFn cfn_add(const CFn* a, const CFn* b, int sign);
CFn cfn_mul(const CFn* a, const CFn* b);
CFn cfn_inv(const CFn* a);                             // a ≠ 0
CFn cfn_derivative(const CFn* a);

// y ∈ Q(i)(x) with y' + f·y = g for a nonzero polynomial f; false when there is none
bool c_solve_rde(const CPoly* f, const CFn* g, CFn* y);

#endif
