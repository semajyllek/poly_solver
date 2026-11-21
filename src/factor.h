#ifndef FACTOR_H
#define FACTOR_H

#include "polynomial.h"
#include "expression.h"

// poly == content · Π factors[i]^multiplicities[i]; factors are primitive
// integer polynomials with positive leading coefficient, irreducible over Q
typedef ps_factorization Factorization;

// square-free decomposition (Yun): poly == lc · Π parts[i]^(i+1), parts monic,
// pairwise coprime, square-free (a part may be the constant 1)
typedef ps_sqfree SquareFreeDecomp;

SquareFreeDecomp square_free_decompose(const Polynomial* poly);
void free_square_free(SquareFreeDecomp* d);

Rational poly_content(const Polynomial* poly);
Polynomial* poly_primitive_part(const Polynomial* poly);
Factorization* factorize(const Polynomial* poly);
void free_factorization(Factorization* f);

char* factorization_to_string(const Factorization* f);

#endif
