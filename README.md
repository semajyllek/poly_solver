# poly_solver

A fast C library and CLI for symbolic polynomial algebra and linear systems with exact rational arithmetic.

## Features

- **Exact arbitrary-precision arithmetic** (GMP): input literals, coefficients and output are never rounded. `1/3` stays `1/3`, `0.1` is `1/10`, and a 30-digit coefficient prints in full.
- **Simplification**: canonicalizes expressions by combining like terms and reducing rational functions.
- **Complete factoring over Q**: Yun square-free decomposition plus Berlekamp–Zassenhaus (Cantor–Zassenhaus mod p, Hensel lifting, recombination), e.g. `x^4 + 4 = (x^2 - 2x + 2)(x^2 + 2x + 2)`.
- **Differentiation & integration**: every rational function over Q can be integrated. Hermite reduction gives the rational part. Linear and quadratic factors get `ln`/`atan` closed forms. Higher-degree factors use **Lazard–Rioboo–Trager**: rational and quadratic log coefficients come out as `ln`/`atan` (Rioboo's continuous arctangents). When the factor's roots are known exactly (binomials `x^n ± c` and rescaled cyclotomics, such as `x^5 - 1` or `x^4 + 2x^2 + 4`), the answer is written with radicals and real `ln`/`atan`, using `cos(kπ/N)` constants where no square-root form exists. Otherwise it is an exact `RootSum`.
- **Equation solving**: exact rational roots, exact radicals for quadratic factors (including complex roots), and certified real-root intervals (Descartes/VCA bisection) for higher degrees.
- **Linear systems**: solve Ax=b, determinant, inverse, rank and RREF with exact rational arithmetic.

## Build

```bash
make        # build test binary
make cli    # build CLI binary
make clean  # remove build artifacts
```

Requires `gcc`, `make` and [GMP](https://gmplib.org/) (`brew install gmp`; the Makefile looks in `/opt/homebrew`).

```bash
make sanitize && ./test_symbolic_solver   # AddressSanitizer + UBSan build of the tests
```

## Usage

```bash
./poly_solver <action> "<expression>"
```

### Polynomial Actions

| Action          | Example                                         | Output                              |
|-----------------|--------------------------------------------------|--------------------------------------|
| `simplify`      | `./poly_solver simplify "5x^3 + 3x^3 - 2x + 4"` | `8x^3 - 2x + 4`                    |
| `factor`        | `./poly_solver factor "x^2 - 5x + 6"`            | `(x - 2)(x - 3)`                   |
| `factor`        | `./poly_solver factor "x^4 + 4"`                 | `(x^2 - 2x + 2)(x^2 + 2x + 2)`     |
| `differentiate` | `./poly_solver differentiate "3x^2 + 4x + 2"`    | `6x + 4`                           |
| `integrate`     | `./poly_solver integrate "6x^2 + 4x"`            | `2x^3 + 2x^2`                      |
| `integrate`     | `./poly_solver integrate "1/((x-1)(x^2+1))"`     | `1/2 * ln(abs(x - 1)) - 1/4 * ln(abs(x^2 + 1)) - 1/2 * atan(x)` |
| `integrate`     | `./poly_solver integrate "x/(x^4+1)"`            | `1/2 * atan(x^2)`                  |
| `integrate`     | `./poly_solver integrate "1/(x^4+1)"`            | `-1/8 * sqrt(2) * ln(abs(x^2 - sqrt(2) * x + 1)) + 1/4 * sqrt(2) * atan(sqrt(2) * x - 1) + 1/8 * sqrt(2) * ln(abs(x^2 + sqrt(2) * x + 1)) + 1/4 * sqrt(2) * atan(sqrt(2) * x + 1)` |
| `solve`         | `./poly_solver solve "x^3 - 6x^2 + 11x - 6"`    | `x = 1, x = 2, x = 3`             |
| `solve`         | `./poly_solver solve "x^2 - x - 1"`              | `x = 1/2 - 1/2*sqrt(5) ~ -0.6180339887, x = 1/2 + 1/2*sqrt(5) ~ 1.6180339887` |
| `solve`         | `./poly_solver solve "x^3 - 2"`                  | `x ~ 1.2599210498, 2 non-real roots of (x^3 - 2)` |

### Linear Algebra Actions

| Action          | Example                                                           | Output                    |
|-----------------|-------------------------------------------------------------------|---------------------------|
| `solve_system`  | `./poly_solver solve_system "[[1,2],[3,4]]" "[5,11]"`             | `x1 = 1, x2 = 2`        |
| `det`           | `./poly_solver det "[[1,2],[3,4]]"`                               | `-2`                      |
| `inverse`       | `./poly_solver inverse "[[1,2],[3,4]]"`                           | `[[-2, 1] [3/2, -1/2]]`  |
| `rank`          | `./poly_solver rank "[[1,2,3],[4,5,6],[7,8,9]]"`                  | `2`                       |
| `rref`          | `./poly_solver rref "[[1,2,3],[4,5,6]]"`                          | `[[1, 0, -1] [0, 1, 2]]` |

Matrix entries can be integers, fractions or decimals of any size (`1/3`, `-5/7`, `0.25`, `1e30`).

### Expression syntax

- Standard operators: `+`, `-`, `*`, `/`, `^`
- Implicit multiplication: `5x^3` is `5 * x^3`
- Unary minus: `-x`, `-(x + 2)`
- Parentheses: `(x + 1)(x - 1)`
- Number literals: integers, decimals and exponents (`12`, `0.75`, `2.5e-3`), all exact
- `-x^2` means `-(x^2)`. Anything left over after a complete expression (e.g. `= 4`) is an error.
- Exponents must be integers. `x^(1/2)` and division by zero are rejected with an error.

### Limits

- When a logarithm's coefficient is an algebraic number of degree 3 or more, the answer is an exact `RootSum(q(t), t -> t * ln(S(t, x)))`: the sum of `t·ln(S(t, x))` over the roots of `q`. This happens only when the roots aren't known in closed form. Example: `1/(x^5+x+3)`, whose roots can't be written with radicals at all (its Galois group is S5).
- Complex roots of irreducible factors of degree 3 or more are counted but not located.

## Tests

```bash
make && ./test_symbolic_solver
```

- 96 unit tests covering rational arithmetic, polynomial algebra, parsing, canonicalization, factoring, solving, calculus and linear systems (including Hilbert matrices).
- Regression tests for each previously found bug: big literals, truncated output, dropped integral terms, lost factor content, `-x^2` parsed as `(-x)^2`, and trailing input silently ignored.
- Seeded property tests that check invariants on generated inputs:
  - `content · Π factors == input` exactly
  - each root interval brackets a sign change
  - the exact derivative of every antiderivative equals the integrand (for `RootSum` terms, compared exactly at rational points via traces over Q(α)); answers with radicals and trig constants are checked numerically in double precision, and to 60 digits by the SymPy cross-check
  - `A·solve(A, A·x) == A·x`
  - parse(print(p)) == p

  Set `POLY_SEED=<n>` to try other seeds.

A differential check against SymPy is also available. It is optional and needs `pip install sympy`:

```bash
make cli && python3 crosscheck_sympy.py --trials 40   # about a minute; --seed N for others
```

It confirms each factorization multiplies back, each root and root count, and each antiderivative by differentiating it. Answers containing `RootSum`, radicals or trig constants are checked at 60 significant digits, including on a separate stream of binomial and cyclotomic denominators.

## Architecture

1. **Rational arithmetic** (`rational_num.c`): GMP `mpq_t`, exact literal parsing, decimal output
2. **Polynomial algebra** (`polynomial.c`): sparse sorted linked list, Euclidean GCD/XGCD, exact Horner evaluation
3. **Modular and integer polynomials** (`poly_zp.c`, `poly_zx.c`): dense F_p[x] and Z[x]
4. **Expression tree** (`expression.c`): tagged union with canonicalization to polynomial/rational form
5. **Parser** (`token.c`, `parser.c`): shunting-yard with implicit multiplication and unary minus
6. **Factoring** (`factor.c`, `zx_factor.c`): Yun square-free decomposition + Berlekamp–Zassenhaus
7. **Roots and solver** (`roots.c`, `solver.c`): Descartes/VCA isolation, exact refinement, radicals
8. **Integration** (`integrate.c`, `partial_frac.c`): Hermite reduction, log/atan terms, partial fractions
9. **Linear algebra** (`matrix.c`): Gaussian elimination over exact rationals

See [INTERNALS.md](INTERNALS.md) for how each algorithm works.

## License

MIT
