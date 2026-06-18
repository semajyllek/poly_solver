# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test

```bash
make              # Build test binary (test_symbolic_solver)
make cli          # Build CLI binary (poly_solver)
make clean        # Remove object files and binaries
./test_symbolic_solver  # Run all tests (assert-based, exits on first failure)
```

CLI usage:
```bash
./poly_solver simplify "5x^3 + 3x^3 - 2x + 4"
./poly_solver factor "x^2 - 5x + 6"
./poly_solver differentiate "3x^2 + 4x + 2"
./poly_solver integrate "6x^2 + 4x"
./poly_solver solve "x^3 - 6x^2 + 11x - 6"
```

## Architecture

Symbolic algebra engine with five layers:

1. **Rational arithmetic** (`rational_num.c/h`) — Exact `Rational` type (`long long` num/den, GCD-normalized). `rat_int_t` typedef enables future GMP swap. Value type passed by value.

2. **Polynomial algebra** (`polynomial.c/h`) — Sparse sorted linked list of `Term` nodes (descending exponent, `Rational` coefficients, `var_id` for future multivariate). Core operations: add, subtract, multiply, divmod (returns quotient + remainder), GCD (Euclidean, non-destructive, monic output), LCM, derivative, integral, Horner evaluation.

3. **Expression tree** (`expression.c/h`) — Tagged union `Expr` with types: `EXPR_CONST`, `EXPR_VAR`, `EXPR_BINOP` (unified binary ops with `BinOp` discriminant), `EXPR_NEG`, `EXPR_POLY`, `EXPR_RATIONAL_FN`. `canonicalize()` converts any expression tree to polynomial or rational function canonical form. `simplify()` = `canonicalize()`.

4. **Parser** (`parser.c/h`) — Recursive descent with implicit multiplication (`5x^3`), unary minus (`-x`), and standard operator precedence. Single-token lookahead via saved input pointer.

5. **Factoring & Solving** (`factor.c/h`, `solver.c/h`) — Rational root theorem + synthetic division for linear factors. Square-free decomposition via `gcd(f, f')`. Quadratic formula for degree-2 factors. Solver dispatches to factorizer, then solves each factor by degree.

## Key Design Invariants

- **Memory ownership**: Caller owns inputs. Every function returns new allocations. No function mutates or frees its inputs.
- **Canonical forms**: All simplification reduces to `EXPR_POLY` or `EXPR_RATIONAL_FN`. Rational functions are GCD-reduced.
- **Single variable for now**: `Term.var_id` field exists for future multivariate but is always 0.

## Code Style
- keep functions small and targeted
- avoid redundant code
- keep function signatures small. more than 3 args should prompt the use of more curried functions or an object to pass the args.
- add one-line comments to code frequently, particularly when complex logic exists
