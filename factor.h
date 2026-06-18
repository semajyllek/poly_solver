#ifndef FACTOR_H
#define FACTOR_H

#include "polynomial.h"
#include "expression.h"

typedef struct {
    Polynomial** factors;
    int* multiplicities;
    int count;
    Rational content;
} Factorization;

Rational poly_content(const Polynomial* poly);
Polynomial* poly_primitive_part(const Polynomial* poly);
Polynomial* synthetic_divide(const Polynomial* poly, const Rational* root);
Factorization* factorize(const Polynomial* poly);
void free_factorization(Factorization* f);

Expr* factorization_to_expr(const Factorization* f);
char* factorization_to_string(const Factorization* f);

#endif
