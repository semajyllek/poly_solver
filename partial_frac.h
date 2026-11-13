#ifndef PARTIAL_FRAC_H
#define PARTIAL_FRAC_H

#include "polynomial.h"
#include "rational_num.h"
#include <stdbool.h>

// One irreducible denominator factor P with multiplicity m contributes
//   numerators[0]/P + numerators[1]/P^2 + ... + numerators[m-1]/P^m,
// each numerator of degree < deg P.
typedef struct {
    Polynomial* factor;       // irreducible over Q, primitive integer coefficients
    int power;                // multiplicity m
    Polynomial** numerators;  // m numerators
    Rational* coeffs;         // linear factors only: the constant numerators (else NULL)
} PFTerm;

typedef struct {
    Polynomial* poly_part; // polynomial part from long division (NULL if proper)
    PFTerm* terms;
    int num_terms;
} PartialFractionDecomp;

PartialFractionDecomp* partial_fractions(const Polynomial* num, const Polynomial* den);
void free_pf_decomp(PartialFractionDecomp* pf);
char* pf_to_string(const PartialFractionDecomp* pf);

#endif
