#ifndef SOLVER_H
#define SOLVER_H

#include "polynomial.h"
#include "rational_num.h"
#include "roots.h"
#include <stdbool.h>

// an irreducible factor of degree >= 2 and what is known exactly about its roots
typedef struct {
    Polynomial* factor;   // irreducible over Q
    int multiplicity;
    RootInterval* real;   // certified isolating intervals of its real roots, ascending
    int num_real;
} IrrationalFactor;

typedef struct {
    Rational* rational_roots;    // exact, repeated by multiplicity
    int num_rational;
    IrrationalFactor* irr;       // one entry per distinct irreducible factor of degree >= 2
    int num_irr;
    // convenience views derived from irr:
    double* irrational_roots;    // approximations of real irrational roots, by multiplicity
    int num_irrational;
    Polynomial** irreducible_factors; // factors that have non-real roots, by multiplicity
    int num_irreducible;
} Solutions;

Solutions* solve_polynomial(const Polynomial* poly);
char* solutions_to_string(const Solutions* sol);
void free_solutions(Solutions* sol);

#endif
