#ifndef SOLVER_H
#define SOLVER_H

#include "polynomial.h"
#include "rational_num.h"
#include <stdbool.h>

typedef struct {
    Rational* rational_roots;
    int num_rational;
    double* irrational_roots;
    int num_irrational;
    Polynomial** irreducible_factors;
    int num_irreducible;
} Solutions;

Solutions* solve_polynomial(const Polynomial* poly);
char* solutions_to_string(const Solutions* sol);
void free_solutions(Solutions* sol);

#endif
