# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test

```bash
make              # Build test binary (test_symbolic_solver)
make cli          # Build CLI binary (poly_solver)
make clean        # Remove object files and binaries
make sanitize     # Rebuild the tests with AddressSanitizer + UBSan
POLY_SEED=7 ./test_symbolic_solver          # property tests with another seed
python3 crosscheck_sympy.py --trials 40     # optional differential test vs SymPy (needs sympy, make cli)
./test_symbolic_solver  # Run all tests (assert-based, exits on first failure)
                        # unit tests: test_symbolic_solver.c; regression + property tests: test_properties.c
```

CLI usage:
```bash
./poly_solver simplify "5x^3 + 3x^3 - 2x + 4"
./poly_solver factor "x^2 - 5x + 6"
./poly_solver differentiate "3x^2 + 4x + 2"
./poly_solver integrate "6x^2 + 4x"
./poly_solver integrate "1/(x^2+1)" 0 1            # definite: exact value ~ decimal
./poly_solver integrate "1/(x^4+1)" -inf inf       # infinite bounds: inf, -inf (or oo, -oo)
./poly_solver solve "x^3 - 6x^2 + 11x - 6"
./poly_solver solve_system "[[1,2],[3,4]]" "[5,11]"
./poly_solver det "[[1,2],[3,4]]"
./poly_solver inverse "[[1,2],[3,4]]"
./poly_solver rank "[[1,2,3],[4,5,6]]"
./poly_solver rref "[[1,2,3],[4,5,6]]"
```

## Architecture

Symbolic algebra engine, layered bottom-up (INTERNALS.md has the details):

1. **Strings and rationals** (`strbuf.c/h`, `rational_num.c/h`): `StrBuf` is the only way output is built; no fixed-size char buffers. `Rational` wraps GMP `mpq_t` (init/clear lifecycle). `rat_parse_literal` parses number text exactly; nothing goes through `double`.

2. **Polynomial algebra** (`polynomial.c/h`): sparse sorted linked list of `Term` nodes (descending exponent, `Rational` coefficients, `var_id` reserved for multivariate). Add, subtract, multiply, divmod, GCD (monic), XGCD, LCM, pow, derivative, integral, exact and double evaluation.

3. **Modular/integer kernels** (`poly_zp.c/h`, `poly_zx.c/h`, `zx_factor.c`): dense F_p[x] (p < 2^31) with distinct-degree and Cantor–Zassenhaus factoring; dense Z[x] with GMP coefficients; Berlekamp–Zassenhaus (prime choice, Hensel lifting, exact-division recombination).

4. **Factoring** (`factor.c/h`): Yun square-free decomposition, then BZ per part. `Factorization` invariant: `poly == content · Π factors[i]^mult[i]`, factors primitive integer polynomials with positive leading coefficient, sorted deterministically.

5. **Roots and solving** (`roots.c/h`, `solver.c/h`): Descartes/VCA real-root isolation with exact bisection refinement. The solver gives exact rational roots, radicals for quadratic factors, every root exactly for binomial/cyclotomic-type factors (`rootform.c`), and certified intervals plus non-real counts for other higher-degree factors.

6. **Expression tree and parser** (`expression.c/h`, `token.c`, `parser.c`): tagged union `Expr` (`CONST`, `VAR`, `BINOP`, `NEG`, `POLY`, `RATIONAL_FN`, `LN`, `ATAN`, `SQRT`, `ROOTSUM`, `TRIG` = cos/sin/cot/csc(kπ/N)). `canonicalize()` reduces to `EXPR_POLY` or `EXPR_RATIONAL_FN` and returns NULL for non-integer exponents or division by zero. Shunting-yard parser with implicit multiplication.

7. **Calculus** (`integrate.c`, `logpart.c/h`, `subres.c/h`, `rootform.c/h`, `surd.c/h`, `partial_frac.c/h`): Hermite reduction (exactly self-checked), closed-form ln/atan for linear and quadratic factors, and Lazard–Rioboo–Trager for higher-degree factors (subresultant PRS chain over Q[t], checked against a determinant oracle in `test_subres_ref.c`; ln, Rioboo LogToAtan, or `EXPR_ROOTSUM`). `integrate()` then rewrites a RootSum into radicals and real ln/atan when the factor is a binomial or a rescaled cyclotomic in x^s (`rootform.c`; exact constants c·m^e·sqrt(s)·trig(kπ/N) in `surd.c`). `integrate_rootsum()` keeps the RootSum form, which is what the exact tests use. Like terms are merged. Partial fractions handle factors of any degree. Definite integrals (`definite.c`): exact pole check, `F(b) − F(a)` folded by `constfold.c`, atan jump corrections, limits at ±∞ from each term's leading behaviour (`expr_lead`), and a run-time comparison with quadrature (numeric-only output when it fails or F has a RootSum).

8. **Linear algebra** (`matrix.c/h`): dense row-major `Matrix` of `Rational`; Gaussian elimination to RREF; solve, det, inverse, rank; bracket-notation parser using the exact literal parser.

## Key Design Invariants

- **Memory ownership**: Caller owns inputs. Every function returns new allocations. No function mutates or frees its inputs.
- **Canonical forms**: All simplification reduces to `EXPR_POLY` or `EXPR_RATIONAL_FN`. Rational functions are GCD-reduced; constant denominators are folded into the polynomial.
- **Exactness**: no `double` on any path that produces a symbolic answer. Doubles appear only in `poly_evaluate`/`expr_eval` (numeric checks), the legacy `Solutions.irrational_roots` view, the decimals printed beside exact complex roots, and `definite.c` (the quadrature check, and reading the side signs of an atan argument next to an exactly isolated pole, or of its leading coefficient at ±∞, both backed up by that check).
- **Never partial**: an operation that cannot finish returns NULL with a reason; it never returns a silently incomplete result.
- **Single variable for now**: `Term.var_id` field exists for future multivariate but is always 0.

## Code Style
- keep functions small and targeted
- avoid redundant code
- keep function signatures small. more than 3 args should prompt the use of more curried functions or an object to pass the args.
- add one-line comments to code frequently, particularly when complex logic exists
