#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <math.h>
#include "rational_num.h"

typedef struct Term {
    Rational coeff;
    int exponent;
    int var_id; // 0 = "x"; future multivariate hook
    struct Term* next;
} Term;

typedef struct Polynomial {
    Term* head;
} Polynomial;

typedef struct {
    Polynomial* quotient;
    Polynomial* remainder;
} PolyDivResult;

// roots of a quadratic: p ± r·sqrt(m) (real) or p ± r·sqrt(m)·i (complex), m square-free-ish
typedef struct {
    Rational p, r;
    mpz_t m;
    bool is_complex;
} QuadRoots;

// s*a + t*b == g with g = gcd(a, b) monic
typedef struct {
    Polynomial* s;
    Polynomial* t;
    Polynomial* g;
} PolyBezout;

// construction
Polynomial* create_polynomial(void);
Term* create_term(const Rational* coeff, int exponent);
void add_term(Polynomial* poly, const Rational* coeff, int exponent);

// copy / query
Polynomial* poly_copy(const Polynomial* poly);
bool poly_is_zero(const Polynomial* poly);
int polynomial_degree(const Polynomial* poly);

// arithmetic
Polynomial* poly_add(const Polynomial* poly1, const Polynomial* poly2);
Polynomial* poly_subtract(const Polynomial* poly1, const Polynomial* poly2);
Polynomial* poly_multiply(const Polynomial* poly1, const Polynomial* poly2);
PolyDivResult poly_divmod(const Polynomial* dividend, const Polynomial* divisor);
Polynomial* poly_divide(const Polynomial* dividend, const Polynomial* divisor);
Polynomial* poly_mod(const Polynomial* dividend, const Polynomial* divisor);
Polynomial* poly_gcd(const Polynomial* poly1, const Polynomial* poly2);
Polynomial* poly_lcm(const Polynomial* poly1, const Polynomial* poly2);
Polynomial* poly_pow(const Polynomial* base, unsigned exp); // binary exponentiation
Polynomial* poly_scale(const Polynomial* p, const Rational* c);
PolyBezout  poly_xgcd(const Polynomial* a, const Polynomial* b);
void        free_bezout(PolyBezout* bz);
Polynomial* poly_constant(const Rational* c);               // the constant polynomial c
void poly_coeff(Rational* out, const Polynomial* poly, int degree); // 0 when absent
void poly_make_monic(Polynomial* poly);

// calculus
Polynomial* poly_derivative(const Polynomial* poly);
Polynomial* poly_integral(const Polynomial* poly);
double poly_definite_integral(const Polynomial* poly, double lower, double upper);

// evaluation and output
double poly_evaluate(const Polynomial* poly, double x);
void poly_eval_rat(Rational* out, const Polynomial* poly, const Rational* x); // exact Horner
char* poly_to_string(const Polynomial* poly); // caller frees
char* poly_to_string_var(const Polynomial* poly, const char* var); // e.g. "t" for resultants
QuadRoots poly_quad_roots(const Polynomial* f); // exact roots of a degree-2 polynomial
void quad_roots_clear(QuadRoots* q);
void free_polynomial(Polynomial* poly);
void print_polynomial(const Polynomial* poly);

#endif
