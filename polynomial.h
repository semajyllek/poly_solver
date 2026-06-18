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

// construction
Polynomial* create_polynomial(void);
Term* create_term(Rational coeff, int exponent);
void add_term(Polynomial* poly, Rational coeff, int exponent);

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
void poly_make_monic(Polynomial* poly);

// calculus
Polynomial* poly_derivative(const Polynomial* poly);
Polynomial* poly_integral(const Polynomial* poly);
double poly_definite_integral(const Polynomial* poly, double lower, double upper);

// evaluation and output
double poly_evaluate(const Polynomial* poly, double x);
char* poly_to_string(const Polynomial* poly); // caller frees
void free_polynomial(Polynomial* poly);
void print_polynomial(const Polynomial* poly);

#endif
