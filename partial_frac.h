#ifndef PARTIAL_FRAC_H
#define PARTIAL_FRAC_H

#include "polynomial.h"
#include "rational_num.h"
#include <stdbool.h>

typedef struct {
    Rational* coeffs;      // A₁..Aₘ for linear, or NULL for quadratic
    Polynomial* factor;    // (x - r) or (x² + bx + c)
    int power;             // highest power m
    bool is_quadratic;
    Rational* quad_coeffs; // B₁,C₁,B₂,C₂,... for quadratic (NULL for linear)
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
