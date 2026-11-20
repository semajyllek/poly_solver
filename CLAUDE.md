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
./poly_solver integrate "x^2*sin(x)" 0 pi           # exp/sin/cos and bounds q·pi: pi^2 - 4
./poly_solver integrate "x^2*exp(3x)"               # Risch: answer, or a proof there is none (exp(x^2))
./poly_solver integrate "x^2*ln(x)^3"               # Risch, ln case (ln(x)/(x+1): proved non-elementary)
./poly_solver integrate "1/(1+exp(x))"              # denominators in θ: Hermite + residue test (1/ln(x): proof)
./poly_solver integrate "x*exp(x)*sin(x)"           # sin/cos as complex exponentials (sin(x)/x: proof)
./poly_solver integrate "1/(2+cos(x))"              # sin/cos denominators: θ = tan(x) or tan(x/2) (x*tan(x): proof)
./poly_solver integrate "exp(x)*(ln(x)+1/x)"        # ln over an exp/ln tower: exp(x)·ln|x| (exp(x)*ln(x): proof)
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

4. **Factoring** (`factor.c/h`, `modgcd.c/h`): Yun square-free decomposition (gcds by `poly_gcd_fast`: CRT over word-size primes + rational reconstruction, checked exactly over Q, Euclid as fallback; also the inverse in LRT's `reduce_monic`), then BZ per part. `Factorization` invariant: `poly == content · Π factors[i]^mult[i]`, factors primitive integer polynomials with positive leading coefficient, sorted deterministically.

5. **Roots and solving** (`roots.c/h`, `solver.c/h`): Descartes/VCA real-root isolation with exact bisection refinement. The solver gives exact rational roots, radicals for quadratic factors, every root exactly for binomial/cyclotomic-type factors (`rootform.c`), and certified intervals plus non-real counts for other higher-degree factors.

6. **Expression tree and parser** (`expression.c/h`, `token.c`, `parser.c`): tagged union `Expr` (`CONST`, `VAR`, `BINOP`, `NEG`, `POLY`, `RATIONAL_FN`, `LN`, `ATAN`, `SQRT`, `ROOTSUM`, `TRIG` = cos/sin/cot/csc(kπ/N), `EXP`, `SIN`, `COS`). Function calls `exp(`, `ln(`, `sin(`, `cos(`, `tan(` (= sin/cos) are tokens; `parse_input` gives the raw tree (`main.c` routes trees with exp/ln/sin/cos to `risch.c`). `canonicalize()` reduces to `EXPR_POLY` or `EXPR_RATIONAL_FN` and returns NULL for non-integer exponents or division by zero. Shunting-yard parser with implicit multiplication.

7. **Calculus** (`integrate.c`, `logpart.c/h`, `subres.c/h`, `rootform.c/h`, `surd.c/h`, `partial_frac.c/h`): Hermite reduction (exactly self-checked), closed-form ln/atan for linear and quadratic factors, and Lazard–Rioboo–Trager for higher-degree factors (subresultant PRS chain over Q[t], checked against a determinant oracle in `test_subres_ref.c`; ln, Rioboo LogToAtan, or `EXPR_ROOTSUM`). `integrate()` then rewrites a RootSum into radicals and real ln/atan when the factor is a binomial or a rescaled cyclotomic in x^s (`rootform.c`; exact constants c·m^e·sqrt(s)·trig(kπ/N) in `surd.c`, square roots merged as sqrt(2)·sqrt(3) = sqrt(6), off-table trig values kept symbolic on purpose). `integrate_rootsum()` keeps the RootSum form, which is what the exact tests use. Like terms are merged. Partial fractions handle factors of any degree. Exp and ln integrands (`expfield.c`, `risch.c`, `kpoly.c`, `ksurd.c`, `residue.c`, `cfield.c`, `expsum.c`, `tanfield.c`, `deftower.c`, `risch2.c`): one tower, θ = exp(v) (v a polynomial; Σ r_k·θ^k, Risch differential equation per k) or θ = ln(u) (u rational; polynomials in θ, limited integration level by level via `hermite_split` in `integrate.h`); a denominator in θ (`exp_form` gives a reduced fraction over K[θ], K = Q(x), arithmetic in `kpoly.c`) goes through Hermite reduction with the derivation and the Rothstein–Trager residue test in `residue.c` (power sums of the residues must be rational; rational residues give c·ln(G_c); quadratic ones a ± b√m give a gcd over K(√M) in `ksurd.c`, then real logs or Rioboo arctangents of polynomials in θ, exactly checked; degree >= 3 is unsupported), and the exact remainder is a polynomial in θ handled as above; sin/cos (and exp-only integrands outside one tower, such as exp(x)·exp(x²)) go to `expsum.c`: exponential sums Σ c_w·e^w over Q(i)(x) (`cfield.c`), w = E + iF polynomials, one complex Risch differential equation per exponent, answers written back with real exp·sin/cos; a sum of sin/cos in a denominator goes to θ = tan(H) (`tanfield.c`: H = B when the integrand is even under B -> B + π, else B/2; derivation θ' = H'(1 + θ²)), where the normal part takes Hermite and the residue test, the polynomial part reduces to p1·θ (p1/(2H') must be a constant: c·ln(1 + θ²) = -2c·ln|cos H|), and the special part over (1 + θ²)^m is a trig polynomial for `expsum.c`; logs of polynomials in tan(H) are rewritten in sin/cos when the ln|cos H| terms cancel; two levels θ2 = ln(u) over K1 = Q(x)(θ1), θ1 = exp or ln, for polynomials in θ2 (`risch2.c`, entered only when the one-tower conversion fails on an integrand with ln): the levels a_k = D(b_k) + (k+1)·b_{k+1}·η with limited integration over K1 by Hermite reduction alone (simple parts in ratio c, θ1^k parts integrated by `exp_part`/`solve_log_levels` exposed in `rischk1.h`, θ1^0 Hermite simple parts s(q_g) = c·s(p_η) + d·L1 solved over Q), each step and the whole checked exactly; each gives the answer or a proof of non-elementarity; anything else is "unsupported", never "non-elementary". Definite integrals (`definite.c`): exact pole check, `F(b) − F(a)` folded by `constfold.c` (ln/atan terms merged and log arguments made canonical with Q(√m) arithmetic in `qsurd.c`), atan jump corrections, limits at ±∞ from each term's leading behaviour (`expr_lead`), and a run-time comparison with quadrature (numeric-only output when it fails or F has a RootSum). Exponential sums (exp/sin/cos, no sum in a denominator; `expsum_shape`) and bounds q·π (`parse_endpoint`) go through `definite_exp`: a coefficient-denominator root in [a, b] is unsupported (sin(x)/x is finite at 0) unless f is rational, an infinite end needs every e^w to decay there (then F's exp/sin/cos terms tend to 0) and c_0 to decay like 1/x², q·π is a 50-digit rational in a band checked free of roots, `expr_at_pi` evaluates exactly (Σ c·π^k, sin/cos(rπ) from the `surd.c` table, quarter-turn shifts sin(a + π) = -sin(a)), and a non-elementary integrand gets the numeric value alone. A sum of exp or sin/cos in a denominator (finite ranges only) goes through `definite_tower`: `deftower.c` finds the poles of N/D in θ = exp(v) or tan(H) (D with constant coefficients, v or H linear; θ = ∞ points and roots r of D, compared with the bounds exactly as a + b·π or in doubles with a margin for atan(r), ln(r)), and the tan(H) crossings inside, where atan(P(tan H)) terms get jump corrections read from the argument's signs and a bound at a crossing takes one-sided limits.

8. **Linear algebra** (`matrix.c/h`): dense row-major `Matrix` of `Rational`; Gaussian elimination to RREF; solve, det, inverse, rank; bracket-notation parser using the exact literal parser.

## Key Design Invariants

- **Memory ownership**: Caller owns inputs. Every function returns new allocations. No function mutates or frees its inputs.
- **Canonical forms**: All simplification reduces to `EXPR_POLY` or `EXPR_RATIONAL_FN`; expressions with `exp` reduce to the Laurent form Σ a_k·exp(k·v) (`simplify_exp`). Rational functions are GCD-reduced; constant denominators are folded into the polynomial.
- **Exactness**: no `double` on any path that produces a symbolic answer. Doubles appear only in `poly_evaluate`/`expr_eval` (numeric checks), the legacy `Solutions.irrational_roots` view, the decimals printed beside exact complex roots, and `definite.c` (the quadrature check, and reading the side signs of an atan argument next to an exactly isolated pole, or of its leading coefficient at ±∞, both backed up by that check).
- **Never partial**: an operation that cannot finish returns NULL with a reason; it never returns a silently incomplete result.
- **Single variable for now**: `Term.var_id` field exists for future multivariate but is always 0.

## Code Style
- keep functions small and targeted
- avoid redundant code
- keep function signatures small. more than 3 args should prompt the use of more curried functions or an object to pass the args.
- add one-line comments to code frequently, particularly when complex logic exists
