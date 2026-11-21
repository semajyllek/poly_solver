#ifndef SOLVER_H
#define SOLVER_H

#include "polynomial.h"
#include "rational_num.h"
#include "roots.h"
#include <stdbool.h>

// an irreducible factor of degree >= 2 and what is known exactly about its roots
typedef ps_irrational_factor IrrationalFactor;

typedef ps_solutions Solutions;

Solutions* solve_polynomial(const Polynomial* poly);
char* solutions_to_string(const Solutions* sol);
void free_solutions(Solutions* sol);

#endif
