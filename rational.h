#ifndef RATIONAL_H
#define RATIONAL_H

#include "polynomial.h" // Include polynomial as it is used

typedef struct RationalExpr {
    Polynomial* numerator;   // Pointer to the numerator polynomial
    Polynomial* denominator; // Pointer to the denominator polynomial
} RationalExpr;

// Function prototypes for rational expressions
RationalExpr* create_rational(Polynomial* numerator, Polynomial* denominator);
void free_rational(RationalExpr* rational);
void print_rational(const RationalExpr* rational);

#endif // RATIONAL_H
