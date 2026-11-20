// polysolver: the whole public API in one include. Link with libpolysolver.a and GMP:
//   make lib && cc -I. my_prog.c libpolysolver.a -lgmp -lm
// Ownership: every function returns new allocations that the caller frees (free_expr,
// free_polynomial, free(char*) for strings, ...), and no function frees or changes its inputs.
#ifndef POLYSOLVER_H
#define POLYSOLVER_H

#include "rational_num.h" // Rational: exact GMP rationals, rat_parse_literal
#include "polynomial.h"   // Polynomial: arithmetic, gcd/xgcd, derivative, exact evaluation
#include "factor.h"       // factorize: complete factoring over Q; square-free decomposition
#include "roots.h"        // isolate_real_roots, refine_root: certified real-root intervals
#include "solver.h"       // solve_polynomial: rational roots, radicals, intervals
#include "expression.h"   // Expr: build, canonicalize, differentiate, integrate (rational), print, eval
#include "parser.h"       // parse_input: text -> Expr
#include "risch.h"        // integrate_exp: exp/ln/sin/cos (Risch), with proofs of non-elementarity
#include "definite.h"     // integrate_definite: exact values with a quadrature check
#include "matrix.h"       // Matrix/Vector: rref, rank, det, inverse, linear systems

#endif
