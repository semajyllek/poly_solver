# poly_solver

A fast C library and CLI for symbolic polynomial algebra with exact rational arithmetic.

## Features

- **Exact rational coefficients** — no floating-point rounding errors; `1/3` stays `1/3`
- **Simplification** — canonicalizes expressions by combining like terms and reducing rational functions
- **Factoring** — rational root theorem, synthetic division, square-free decomposition
- **Differentiation & Integration** — symbolic calculus on polynomial expressions
- **Equation solving** — finds rational roots exactly, irrational roots numerically

## Build

```bash
make        # build test binary
make cli    # build CLI binary
make clean  # remove build artifacts
```

Requires `gcc` and `make`. No external dependencies.

## Usage

```bash
./poly_solver <action> "<expression>"
```

### Actions

| Action          | Example                                         | Output                              |
|-----------------|--------------------------------------------------|--------------------------------------|
| `simplify`      | `./poly_solver simplify "5x^3 + 3x^3 - 2x + 4"` | `8x^3 - 2x + 4`                    |
| `factor`        | `./poly_solver factor "x^2 - 5x + 6"`            | `(x - 2)(x - 3)`                   |
| `differentiate` | `./poly_solver differentiate "3x^2 + 4x + 2"`    | `6x + 4`                           |
| `integrate`     | `./poly_solver integrate "6x^2 + 4x"`            | `2x^3 + 2x^2`                      |
| `solve`         | `./poly_solver solve "x^3 - 6x^2 + 11x - 6"`    | `x = 1, x = 2, x = 3`             |

### Expression syntax

- Standard operators: `+`, `-`, `*`, `/`, `^`
- Implicit multiplication: `5x^3` is `5 * x^3`
- Unary minus: `-x`, `-(x + 2)`
- Parentheses: `(x + 1)(x - 1)`

## Tests

```bash
make && ./test_symbolic_solver
```

42 tests covering rational arithmetic, polynomial algebra, parsing, canonicalization, factoring, solving, and calculus.

## Architecture

1. **Rational arithmetic** (`rational_num.c`) — exact `num/den` with GCD normalization
2. **Polynomial algebra** (`polynomial.c`) — sparse sorted linked list, Euclidean GCD, Horner evaluation
3. **Expression tree** (`expression.c`) — tagged union with canonicalization to polynomial/rational form
4. **Parser** (`parser.c`) — recursive descent with implicit multiplication and unary minus
5. **Factoring** (`factor.c`) — rational root theorem + synthetic division
6. **Solver** (`solver.c`) — dispatches to factorizer, solves each factor by degree

## License

MIT
