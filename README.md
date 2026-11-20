# poly_solver

A fast C library and CLI for exact symbolic algebra: simplification, factoring, equation solving, differentiation, indefinite and definite integration (rational functions, and exp, ln, sin and cos through the Risch algorithm), and linear algebra. Every number is an arbitrary-precision GMP rational, so nothing is ever rounded. When an integral has no elementary antiderivative, you get a proof of that, not a failure.

```console
$ ./poly_solver factor "x^4 + 4"
(x^2 - 2x + 2)(x^2 + 2x + 2)

$ ./poly_solver integrate "x^2*exp(3x)"
(1/3x^2 - 2/9x + 2/27) * exp(3x)

$ ./poly_solver integrate "exp(x^2)"
no elementary antiderivative: the part exp(x^2) (Risch: y' + k·v'·y = a has no rational solution)
```

## Contents

- [Features](#features)
- [Build](#build)
- [Command-line cookbook](#command-line-cookbook)
- [Using the C library](#using-the-c-library)
- [Expression syntax](#expression-syntax)
- [Limits](#limits)
- [Tests](#tests)
- [Architecture](#architecture)

## Features

- **Exact arithmetic of any size** (GMP). Literals are parsed exactly (`0.1` is `1/10`, `2.5e30` keeps every digit), and coefficients, roots, determinants and integrals never overflow or round.
- **Simplification** to a canonical polynomial or reduced rational function, and a canonical form for exp, sin and cos (`sin(x)^2 + cos(x)^2` is `1`).
- **Complete factoring over Q**: Yun square-free decomposition, then Berlekamp–Zassenhaus (Cantor–Zassenhaus mod p, Hensel lifting, recombination).
- **Equation solving**: exact rational roots and radicals for quadratics, including complex roots. Binomial and cyclotomic-type factors (`x^6 + 2`, `x^5 - 1`) get every root exactly. Other factors get certified real-root intervals (Descartes/VCA) and a count of their non-real roots.
- **Rational integration**: Hermite reduction plus Lazard–Rioboo–Trager. The logs and arctangents come out real. Radicals are used when the roots are known, and an exact `RootSum` otherwise.
- **Risch integration** with exp, ln, sin, cos and tan. It covers one exponential or logarithmic tower, sums of exponentials, sin/cos in a denominator (through `tan(x)` or `tan(x/2)`), and ln over an exp or ln tower. Every answer is checked exactly, and "no elementary antiderivative" is always a proof.
- **Definite integrals**:
  - exact values with constants folded (`ln(4/3)`, `1/9 * sqrt(3) * pi`);
  - infinite bounds and bounds at rational multiples of π;
  - divergence detected exactly, and jumps of arctangents corrected;
  - every exact value compared with numeric quadrature at run time.
- **Linear algebra over Q**: solve `Ax = b` (unique, none or infinitely many solutions), determinant, inverse, rank and RREF.

## Build

```bash
make              # the test binary (test_symbolic_solver)
make cli          # the command-line tool (poly_solver)
make lib          # the static library (libpolysolver.a)
make examples     # build and run the C examples in examples/
make check        # the suite plus the golden CLI corpus (a few seconds)
make clean
```

You need `gcc` (or `clang`), `make` and [GMP](https://gmplib.org/): `brew install gmp` on macOS, or `apt install libgmp-dev` on Debian/Ubuntu. The Makefile looks in `/opt/homebrew`; adjust `CFLAGS`/`LDFLAGS` for other prefixes.

## Command-line cookbook

```bash
./poly_solver <action> "<expression>"
./poly_solver integrate "<expression>" <a> <b>        # definite: bounds are rationals, q·pi, inf or -inf
./poly_solver solve_system "<matrix>" "<vector>"      # also det, inverse, rank, rref "<matrix>"
```

Every output below was produced by the current build. A trailing `~ 1.234` is a decimal beside an exact value. Answers write `ln(abs(u))`, which is valid on both sides of a sign change of `u`.

### Simplify

Like terms combine, rational functions are reduced by their gcd, decimals become exact fractions, and exponentials and sin/cos reach a canonical form.

```console
$ ./poly_solver simplify "5x^3 + 3x^3 - 2x + 4"
8x^3 - 2x + 4

$ ./poly_solver simplify "(x+1)^3 - (x-1)^3"
6x^2 + 2

$ ./poly_solver simplify "(x^2-1)/(x^2+2x+1)"
(x - 1) / (x + 1)

$ ./poly_solver simplify "0.25x + 1/3x"
7/12x

$ ./poly_solver simplify "(exp(2x) - 1)/(exp(x) - 1)"
1 + exp(x)

$ ./poly_solver simplify "sin(x)^2 + cos(x)^2"
1
```

### Factor

The factorization is complete over Q. The content stays in front, and factors are primitive integer polynomials sorted deterministically.

```console
$ ./poly_solver factor "x^2 - 5x + 6"
(x - 2)(x - 3)

$ ./poly_solver factor "6x^4 + 5x^3 - 38x^2 + 5x + 6"
(x + 3)(3x + 1)(2x - 1)(x - 2)

$ ./poly_solver factor "x^4 + 4"
(x^2 - 2x + 2)(x^2 + 2x + 2)

$ ./poly_solver factor "x^6 - 1"
(x + 1)(x - 1)(x^2 - x + 1)(x^2 + x + 1)

$ ./poly_solver factor "x^8 + 4"
(x^4 - 2x^2 + 2)(x^4 + 2x^2 + 2)

$ ./poly_solver factor "2x^6 - 2x^4 + 8x^2 - 8"
2(x + 1)(x - 1)(x^2 - 2x + 2)(x^2 + 2x + 2)

$ ./poly_solver factor "(x^2+1)^20*(x-3)^7"
(x - 3)^7(x^2 + 1)^20
```

### Solve

The equation is `p(x) = 0`. Rational roots are exact. Quadratic factors give radicals, real or complex. Binomials and cyclotomic-type factors give every root exactly, with `cos(kπ/N)` where no square root form exists. Other factors give certified decimals (each interval is proved to hold exactly one root) and a count of their non-real roots.

```console
$ ./poly_solver solve "x^3 - 6x^2 + 11x - 6"
x = 1, x = 2, x = 3

$ ./poly_solver solve "x^2 - x - 1"
x = 1/2 - 1/2*sqrt(5) ~ -0.6180339887, x = 1/2 + 1/2*sqrt(5) ~ 1.6180339887

$ ./poly_solver solve "x^2 - 2x + 5"
x = 1 +/- 2*i

$ ./poly_solver solve "x^3 - 2"
x = 2^(1/3) ~ 1.2599210498, x = -1/2 * 2^(1/3) +/- 1/2 * 2^(1/3) * sqrt(3) * i ~ (-0.6299605249 +/- 1.0911236360i)

$ ./poly_solver solve "x^5 - 1"
x = 1, x = cos(2*pi/5) +/- sin(2*pi/5) * i ~ (0.3090169944 +/- 0.9510565163i), x = -cos(pi/5) +/- sin(pi/5) * i ~ (-0.8090169944 +/- 0.5877852523i)

$ ./poly_solver solve "x^6 + 2"
x = 1/2 * 2^(1/6) * sqrt(3) +/- 1/2 * 2^(1/6) * i ~ (0.9720806486 +/- 0.5612310242i), x = +/- 2^(1/6) * i ~ (0.0000000000 +/- 1.1224620483i), x = -1/2 * 2^(1/6) * sqrt(3) +/- 1/2 * 2^(1/6) * i ~ (-0.9720806486 +/- 0.5612310242i)

$ ./poly_solver solve "x^4 - 10x^2 + 1"
x ~ -3.1462643699, x ~ -0.3178372451, x ~ 0.3178372451, x ~ 3.1462643699

$ ./poly_solver solve "x^3 + x + 1"
x ~ -0.6823278038, 2 non-real roots of (x^3 + x + 1)
```

### Differentiate

```console
$ ./poly_solver differentiate "3x^2 + 4x + 2"
6x + 4

$ ./poly_solver differentiate "x^3/(x^2+1)"
(x^4 + 3x^2) / (x^4 + 2x^2 + 1)

$ ./poly_solver differentiate "x^2*exp(3x)"
(3x^2 + 2x) * exp(3x)

$ ./poly_solver differentiate "exp(x^2)*sin(3x)"
2x * exp(x^2) * sin(3x) + 3 * exp(x^2) * cos(3x)

$ ./poly_solver differentiate "ln(x^2+1)^2"
4x * ln(abs(x^2 + 1)) / (x^2 + 1)

$ ./poly_solver differentiate "tan(x)"
tan(x)^2 + 1
```

### Integrate rational functions

Every rational function integrates. The answer uses logs and arctangents with real arguments. Radicals appear when the factor's roots are known in closed form, and a `RootSum` otherwise. `RootSum(q(t), t -> t * ln(S))` is the sum of `t·ln(S)` over the roots `t` of `q`; `x^5 + x + 3` has no roots in radicals at all (its Galois group is S5).

```console
$ ./poly_solver integrate "6x^2 + 4x"
2x^3 + 2x^2

$ ./poly_solver integrate "(3x^2+1)/(x^3+x)"
ln(abs(x)) + ln(abs(x^2 + 1))

$ ./poly_solver integrate "1/((x-1)(x^2+1))"
1/2 * ln(abs(x - 1)) - 1/4 * ln(abs(x^2 + 1)) - 1/2 * atan(x)

$ ./poly_solver integrate "x/(x^4+1)"
1/2 * atan(x^2)

$ ./poly_solver integrate "1/(x^2-2)"
1/4 * sqrt(2) * ln(abs((x - sqrt(2)) / (x + sqrt(2))))

$ ./poly_solver integrate "1/(x^4+1)"
-1/8 * sqrt(2) * ln(abs(x^2 - sqrt(2) * x + 1)) + 1/4 * sqrt(2) * atan(sqrt(2) * x - 1) + 1/8 * sqrt(2) * ln(abs(x^2 + sqrt(2) * x + 1)) + 1/4 * sqrt(2) * atan(sqrt(2) * x + 1)

$ ./poly_solver integrate "1/(x^3-2)"
1/6 * 2^(1/3) * ln(abs(x - 2^(1/3))) - 1/12 * 2^(1/3) * ln(abs(x^2 + 2^(1/3) * x + 2^(2/3))) - 1/6 * 2^(1/3) * sqrt(3) * atan(1/3 * 2^(2/3) * sqrt(3) * x + 1/3 * sqrt(3))

$ ./poly_solver integrate "1/(x^5+x+3)"
RootSum(253381t^5 - 160t^3 - 80t^2 - 15t - 1, t -> t * ln(x + (-16216384/16875t^4 + 4054096/16875t^3 - 111476/1875t^2 + 255941/16875t + 256/16875)))
```

### Integrate with exp

Here `r(x)·exp(k·v)` sums are handled, `v` a polynomial, one Risch differential equation per power. Products of different exponentials merge. A "no elementary antiderivative" result is a proof, and it names the part that fails.

```console
$ ./poly_solver integrate "x^2*exp(3x)"
(1/3x^2 - 2/9x + 2/27) * exp(3x)

$ ./poly_solver integrate "x^3*exp(-x^2)"
(-1/2x^2 - 1/2) * exp(-x^2)

$ ./poly_solver integrate "(x-1)*exp(x)/x^2"
exp(x) / x

$ ./poly_solver integrate "exp(x)*exp(x^2)*(2x+1)"
exp(x^2 + x)

$ ./poly_solver integrate "exp(x)/x"
no elementary antiderivative: the part exp(x) / x (Risch: y' + k·v'·y = a has no rational solution)

$ ./poly_solver integrate "exp(x^2)"
no elementary antiderivative: the part exp(x^2) (Risch: y' + k·v'·y = a has no rational solution)
```

### Integrate with ln

```console
$ ./poly_solver integrate "ln(x)^2"
x * ln(abs(x))^2 - 2x * ln(abs(x)) + 2x

$ ./poly_solver integrate "x^2*ln(x)^3"
1/3x^3 * ln(abs(x))^3 - 1/3x^3 * ln(abs(x))^2 + 2/9x^3 * ln(abs(x)) - 2/27x^3

$ ./poly_solver integrate "ln(x)/x^2"
-ln(abs(x)) / x - 1 / x

$ ./poly_solver integrate "ln(x^2+1)"
x * ln(abs(x^2 + 1)) - 2x + 2 * atan(x)

$ ./poly_solver integrate "ln(x)/(x+1)"
no elementary antiderivative: the part ln(abs(x)) / (x + 1) (Risch: a ln level needs a logarithm other than the one present)
```

### exp or ln in a denominator

These go through Hermite reduction and the residue test. Rational residues give logs, and quadratic ones give logs or arctangents of polynomials in the tower. Residues of degree 3 to 6 give an exact RootSum over the tower. Residues that aren't constant are a proof of non-elementarity (`1/ln(x)` is the logarithmic integral).

```console
$ ./poly_solver integrate "1/(1+exp(x))"
x - ln(exp(x) + 1)

$ ./poly_solver integrate "exp(x)/(1+exp(2x))"
atan(exp(x))

$ ./poly_solver integrate "1/(x*ln(x))"
ln(abs(ln(abs(x))))

$ ./poly_solver integrate "1/(x*(ln(x)^2+3))"
1/3 * sqrt(3) * atan(1/3 * sqrt(3) * ln(abs(x)))

$ ./poly_solver integrate "exp(x)/(exp(3x)+2)"
RootSum(108t^3 - 1, t -> t * ln(exp(x) + 6 * t))

$ ./poly_solver integrate "1/ln(x)"
no elementary antiderivative: the part 1 / ln(abs(x)) (Risch: its residues are not constants)

$ ./poly_solver integrate "x/(1+exp(x))"
no elementary antiderivative: the part x / (exp(x) + 1) (Risch: its residues are not constants)
```

### sin, cos and tan

Sums of exponentials cover the numerators, and a sin/cos denominator is rewritten in θ = `tan(x)` when Bioche's rule allows it, else `tan(x/2)`.

```console
$ ./poly_solver integrate "x*sin(x)"
sin(x) - x * cos(x)

$ ./poly_solver integrate "exp(x)*sin(2x)"
1/5 * exp(x) * sin(2x) - 2/5 * exp(x) * cos(2x)

$ ./poly_solver integrate "sin(x)^3*cos(x)^2"
-1/8 * cos(x) - 1/48 * cos(3x) + 1/80 * cos(5x)

$ ./poly_solver integrate "tan(x)"
-ln(abs(cos(x)))

$ ./poly_solver integrate "1/(2+cos(x))"
2/3 * sqrt(3) * atan(1/3 * sqrt(3) * tan(1/2x))

$ ./poly_solver integrate "1/(1+sin(x))"
-2 / (tan(1/2x) + 1)

$ ./poly_solver integrate "1/(sin(x)+cos(x)^2)"
RootSum(25t^4 - 5t^2 - 1, t -> t * ln(tan(1/2x) + t + 5 * t^2 - 10 * t^3))

$ ./poly_solver integrate "sin(x)/x"
no elementary antiderivative: the part sin(x) / x (Risch: y' + w'·y = c has no rational solution over Q(i))

$ ./poly_solver integrate "x*tan(x)"
no elementary antiderivative: the part x * tan(x) (Risch: p·tan(H) needs p/(2H') constant)
```

### ln over an exp or ln tower

These are two levels, rational in the outer ln. This works for `integrate` only, not `differentiate` or `simplify`.

```console
$ ./poly_solver integrate "exp(x)*(ln(x)+1/x)"
exp(x) * ln(abs(x))

$ ./poly_solver integrate "ln(ln(x))/x"
ln(abs(x)) * ln(abs(ln(abs(x)))) - ln(abs(x))

$ ./poly_solver integrate "ln(x)+ln(x+1)"
x * ln(abs(x + 1)) + x * ln(abs(x)) - 2x + ln(abs(x + 1))

$ ./poly_solver integrate "1/(x*ln(x)*ln(ln(x)))"
ln(abs(ln(abs(ln(abs(x))))))

$ ./poly_solver integrate "1/(x*ln(x)*(ln(ln(x))^2+1))"
atan(ln(abs(ln(abs(x)))))

$ ./poly_solver integrate "1/(x*ln(x)*(ln(ln(x))^3-2))"
RootSum(108t^3 - 1, t -> t * ln(ln(abs(ln(abs(x)))) - 6 * t))

$ ./poly_solver integrate "exp(x)*ln(x)"
no elementary antiderivative: the part -exp(x) / x (Risch: y' + k·v'·y = a has no rational solution)

$ ./poly_solver integrate "ln(ln(x))"
no elementary antiderivative: the part -1 / ln(abs(x)) (Risch: its residues are not constants)
```

### Definite integrals of rational functions

The result is exact, with logs merged and arctangents combined in Q(√m), plus a decimal checked against adaptive quadrature. Poles in range and too-slow decay at infinity are reported.

```console
$ ./poly_solver integrate "1/(x^2+1)" 0 1
1/4 * pi ~ 0.7853981634

$ ./poly_solver integrate "1/(x^2+x+1)" 0 1
1/9 * sqrt(3) * pi ~ 0.6045997881

$ ./poly_solver integrate "1/(x^4+1)" 0 1
1/4 * sqrt(2) * ln(1 + sqrt(2)) + 1/8 * sqrt(2) * pi ~ 0.8669729873

$ ./poly_solver integrate "1/(x^4+1)" -inf inf
1/2 * sqrt(2) * pi ~ 2.221441469

$ ./poly_solver integrate "x/(x^2+1)^2" 0 inf
1/2 ~ 0.5

$ ./poly_solver integrate "1/x" -1 1
Error: the integrand has a pole in [a, b], so the integral diverges.

$ ./poly_solver integrate "1/x" 1 inf
Error: the integral diverges at infinity (the denominator's degree must exceed the numerator's by 2 or more).
```

### Definite integrals with exp, sin and cos

Bounds may be rational multiples of π (`pi`, `pi/2`, `3pi/4`, `2pi`, `-pi`). sin and cos of rational multiples of π take exact values, and other constants stay symbolic (`exp(1)`). For `1/(2+cos(x))` over a full period, the jump of `atan(tan(x/2)/√3)` at π is corrected. Without an elementary antiderivative, the numeric value is given alone.

```console
$ ./poly_solver integrate "x^2*sin(x)" 0 pi
pi^2 - 4 ~ 5.869604401

$ ./poly_solver integrate "sin(x)^2" 0 pi
1/2 * pi ~ 1.570796327

$ ./poly_solver integrate "exp(-x)*cos(x)" 0 inf
1/2 ~ 0.5

$ ./poly_solver integrate "x*exp(-x)" 0 inf
1 ~ 1

$ ./poly_solver integrate "exp(x)" 0 1
exp(1) - 1 ~ 1.718281828

$ ./poly_solver integrate "1/(2+cos(x))" 0 2pi
2/3 * sqrt(3) * pi ~ 3.627598728

$ ./poly_solver integrate "1/(1+exp(x))" 0 1
-ln(1 + exp(1)) + ln(2) + 1 ~ 0.379885493

$ ./poly_solver integrate "1/(1+cos(x))" 0 pi
Error: the integrand has a pole in [a, b], so the integral diverges.

$ ./poly_solver integrate "exp(-x^2)" -inf inf
~ 1.772453851 (no elementary antiderivative, so only the numeric value is given)
```

### Linear algebra

Matrices use bracket notation. Entries may be integers, fractions or decimals of any size, and all arithmetic is exact.

```console
$ ./poly_solver solve_system "[[1,2],[3,4]]" "[5,11]"
x1 = 1, x2 = 2

$ ./poly_solver solve_system "[[1,1,1],[1,2,3],[1,4,9]]" "[6,14,36]"
x1 = 1, x2 = 2, x3 = 3

$ ./poly_solver solve_system "[[1,2],[2,4]]" "[3,6]"
Infinitely many solutions (rank = 1)

$ ./poly_solver solve_system "[[1,2],[2,4]]" "[3,7]"
No solution (inconsistent system, rank = 1)

$ ./poly_solver det "[[1,2],[3,4]]"
-2

$ ./poly_solver det "[[1/2, 1/3], [1/4, 1/5]]"
1/60

$ ./poly_solver inverse "[[2, 1, 1], [1, 3, 2], [1, 0, 0]]"
[[0, 0, 1]
 [-2, 1, 3]
 [3, -1, -5]]

$ ./poly_solver inverse "[[1,2],[2,4]]"
Matrix is singular (no inverse).

$ ./poly_solver rank "[[1,2,3],[4,5,6],[7,8,9]]"
2

$ ./poly_solver rref "[[2,4,-2,2],[4,9,-3,8],[-2,-3,7,10]]"
[[1, 0, 0, -1]
 [0, 1, 0, 2]
 [0, 0, 1, 2]]
```

### Very large (and very small) numbers

Nothing is limited by machine precision: 20- and 40-digit numbers factor and solve exactly, a determinant of `1e300` against `1e-300` is exactly 0, and `2` vs `4.0000000000000000000000001` is a genuine rank difference. Roots near zero are refined to relative precision and printed in scientific form.

```console
$ ./poly_solver simplify "(2x + 0.000000000000000000001)^2"
4x^2 + 1/250000000000000000000x + 1/1000000000000000000000000000000000000000000

$ ./poly_solver factor "(x - 98765432109876543210)^2*(x+1)"
(x + 1)(x - 98765432109876543210)^2

$ ./poly_solver factor "x^128 - 1"
(x + 1)(x - 1)(x^2 + 1)(x^4 + 1)(x^8 + 1)(x^16 + 1)(x^32 + 1)(x^64 + 1)

$ ./poly_solver solve "(x - 12345678901234567890/7)(x + 1)"
x = -1, x = 12345678901234567890/7

$ ./poly_solver solve "x^2 - 100000000000000000000000000000000000000001"
x = -sqrt(100000000000000000000000000000000000000001) ~ -316227766016837933199.8893544432, x = sqrt(100000000000000000000000000000000000000001) ~ 316227766016837933199.8893544432

$ ./poly_solver solve "x^3 - 1000000000000000000000000000000x + 1"
x ~ -1000000000000000.0000000000, x ~ 1.000000000e-30, x ~ 999999999999999.9999999999

$ ./poly_solver integrate "123456789012345678901234567890*x^2*exp(98765x)"
(24691357802469135780246913578/19753x^2 - 49382715604938271560493827156/1950905045x + 49382715604938271560493827156/192681136769425) * exp(98765x)

$ ./poly_solver integrate "x^99999" 0 1
1/100000 ~ 1e-05

$ ./poly_solver integrate "1/(x^2+1)" 0 "1e40"
atan(10000000000000000000000000000000000000000) ~ 1.570796327

$ ./poly_solver integrate "1/(x^2+1e-300)" 0 1
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 * atan(1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000) ~ 1.570796327e+150

$ ./poly_solver integrate "x^2" 0 "1e200"
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000/3 (beyond double precision, so not compared with quadrature)

$ ./poly_solver det "[[1e300, 1], [1, 1e-300]]"
0

$ ./poly_solver det "[[123456789012345678901, 2], [3, 4]]"
493827156049382715598

$ ./poly_solver rank "[[1,2],[2,4.0000000000000000000000001]]"
2
```

## Using the C library

`make lib` builds `libpolysolver.a`. Include the umbrella header `poly_solver.h` (it pulls in every public header) and link with GMP:

```bash
make lib
cc -I/path/to/poly_solver my_prog.c /path/to/poly_solver/libpolysolver.a -lgmp -lm
```

Conventions that hold everywhere:

- **Ownership.** The caller owns its inputs, and every function returns new allocations. Nothing frees or changes its arguments. Free results with `free_expr`, `free_polynomial`, `free_factorization`, `free_solutions`, `matrix_free`, `definite_free`, and so on; strings from `*_to_string` with `free`.
- **Rationals.** `Rational` wraps `mpq_t`. Call `rat_init` before use and `rat_clear` after. The arithmetic functions write into an initialized result (`rat_add(&r, &a, &b)`), and `r->val` is the raw `mpq_t` when you need GMP directly.
- **Never partial.** An operation that cannot finish returns `NULL` (or an error status) with a reason. It never returns a silently incomplete result.
- **One variable.** Expressions are in `x`.

Each program below is a complete file in [`examples/`](examples/); `make examples` builds and runs them all, and the output shown is theirs.

### Exact rationals

```c
// examples/rationals.c
// Exact rationals: literals are parsed exactly (0.1 is 1/10), and nothing is rounded.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, const Rational* r) {
    char* s = rat_to_string(r);
    printf("%-22s %s\n", label, s);
    free(s);
}

int main(void) {
    Rational a, b, sum, big;
    rat_init(&a); rat_init(&b); rat_init(&sum); rat_init(&big);
    rat_parse_literal(&a, "0.1");                  // exactly 1/10
    rat_set_si(&b, 1, 3);
    rat_add(&sum, &a, &b);
    show("0.1 + 1/3 =", &sum);                     // 13/30
    rat_parse_literal(&big, "2.5e30");
    rat_mul(&big, &big, &big);
    show("(2.5e30)^2 =", &big);                    // every digit kept
    char* dec = rat_to_decimal(&b, 20);
    printf("%-22s %s\n", "1/3 to 20 places:", dec);
    free(dec);
    rat_clear(&a); rat_clear(&b); rat_clear(&sum); rat_clear(&big);
    return 0;
}
```

Output:

```text
0.1 + 1/3 =            13/30
(2.5e30)^2 =           6250000000000000000000000000000000000000000000000000000000000
1/3 to 20 places:      0.33333333333333333333
```

### Polynomials: arithmetic, gcd, Bezout, calculus

`Polynomial` is a sparse list of terms with `Rational` coefficients. Build one from text with `parse_input` plus `canonicalize`, or term by term with `add_term`.

```c
// examples/polynomials.c
// Polynomials over Q: build from text or term by term, then arithmetic, gcd and calculus.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

// text -> polynomial (NULL if the text is not a polynomial)
static Polynomial* poly_of(const char* text) {
    Expr* raw = parse_input(text);
    Expr* e = raw ? canonicalize(raw) : NULL;
    Polynomial* p = e && e->type == EXPR_POLY ? poly_copy(e->poly) : NULL;
    free_expr(raw); free_expr(e);
    return p;
}

static void show(const char* label, const Polynomial* p) {
    char* s = poly_to_string(p);
    printf("%-16s %s\n", label, s);
    free(s);
}

int main(void) {
    Polynomial* a = poly_of("(x - 1)^2 * (x + 2)");
    Polynomial* b = create_polynomial();          // x^2 - 1, term by term
    Rational one = rat_from_int(1), minus_one = rat_from_int(-1);
    add_term(b, &one, 2);
    add_term(b, &minus_one, 0);

    show("a =", a);
    show("b =", b);
    Polynomial* prod = poly_multiply(a, b);
    show("a * b =", prod);
    PolyDivResult qr = poly_divmod(a, b);
    show("a div b =", qr.quotient);
    show("a mod b =", qr.remainder);
    Polynomial* g = poly_gcd(a, b);
    show("gcd(a, b) =", g);
    PolyBezout bz = poly_xgcd(a, b);               // s·a + t·b = gcd
    show("Bezout s =", bz.s);
    show("Bezout t =", bz.t);
    Polynomial* da = poly_derivative(a);
    show("a' =", da);

    Rational x0, v;                                // exact evaluation at x = 1/2
    rat_init(&x0); rat_init(&v);
    rat_set_si(&x0, 1, 2);
    poly_eval_rat(&v, a, &x0);
    char* vs = rat_to_string(&v);
    printf("%-16s %s\n", "a(1/2) =", vs);
    free(vs);

    rat_clear(&one); rat_clear(&minus_one); rat_clear(&x0); rat_clear(&v);
    free_bezout(&bz);
    free_polynomial(a); free_polynomial(b); free_polynomial(prod); free_polynomial(g); free_polynomial(da);
    free_polynomial(qr.quotient); free_polynomial(qr.remainder);
    return 0;
}
```

Output:

```text
a =              x^3 - 3x + 2
b =              x^2 - 1
a * b =          x^5 - 4x^3 + 2x^2 + 3x - 2
a div b =        x
a mod b =        -2x + 2
gcd(a, b) =      x - 1
Bezout s =       -1/2
Bezout t =       1/2x
a' =             3x^2 - 3
a(1/2) =         5/8
```

### Factoring, solving and certified roots

`factorize` gives `poly == content · Π factors[i]^multiplicities[i]`. `solve_polynomial` collects the exact and certified roots. `isolate_real_roots` and `refine_root` give intervals with exact rational endpoints, each proved to hold exactly one root, refined as far as you ask.

```c
// examples/factor_solve.c
// Factoring over Q, and solving: exact roots where they exist, certified intervals elsewhere.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static Polynomial* poly_of(const char* text) {
    Expr* raw = parse_input(text);
    Expr* e = raw ? canonicalize(raw) : NULL;
    Polynomial* p = e && e->type == EXPR_POLY ? poly_copy(e->poly) : NULL;
    free_expr(raw); free_expr(e);
    return p;
}

int main(void) {
    Polynomial* p = poly_of("2x^6 - 2x^4 + 8x^2 - 8");
    Factorization* f = factorize(p);               // p == content · Π factors^mult
    char* fs = factorization_to_string(f);
    printf("factor: %s\n", fs);
    free(fs);
    for (int i = 0; i < f->count; i++) {
        char* s = poly_to_string(f->factors[i]);
        printf("  factor %d: %s (multiplicity %d)\n", i + 1, s, f->multiplicities[i]);
        free(s);
    }

    Solutions* sol = solve_polynomial(p);
    char* ss = solutions_to_string(sol);
    printf("solve:  %s\n", ss);
    free(ss);

    Polynomial* q = poly_of("x^5 - 3x + 1");     // no radicals: certified real-root intervals
    RootIntervalList L = isolate_real_roots(q);
    Rational w;
    rat_init(&w);
    rat_set_si(&w, 1, 1000000);
    for (int i = 0; i < L.count; i++) {
        refine_root(&L.items[i], q, &w);          // shrink to width 1e-6, exactly
        char *lo = rat_to_decimal(&L.items[i].lo, 8), *hi = rat_to_decimal(&L.items[i].hi, 8);
        printf("  root of x^5 - 3x + 1 in (%s, %s)\n", lo, hi);
        free(lo); free(hi);
    }
    rat_clear(&w);
    root_interval_list_free(&L);
    free_solutions(sol); free_factorization(f);
    free_polynomial(p); free_polynomial(q);
    return 0;
}
```

Output:

```text
factor: 2(x + 1)(x - 1)(x^2 - 2x + 2)(x^2 + 2x + 2)
  factor 1: x + 1 (multiplicity 1)
  factor 2: x - 1 (multiplicity 1)
  factor 3: x^2 - 2x + 2 (multiplicity 1)
  factor 4: x^2 + 2x + 2 (multiplicity 1)
solve:  x = -1, x = 1, x = 1 +/- i, x = -1 +/- i
  root of x^5 - 3x + 1 in (-1.38879203, -1.38879108)
  root of x^5 - 3x + 1 in (0.33473396, 0.33473491)
  root of x^5 - 3x + 1 in (1.21464729, 1.21464824)
```

### Building and evaluating expressions

Constructors (`create_add`, `create_pow`, `create_var`, ...) take ownership of their arguments. `simplify`, `differentiate` and `integrate` do not; they return new trees.

```c
// examples/expressions.c
// Building expressions in code instead of parsing text, then simplifying, differentiating
// and evaluating them.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, const Expr* e) {
    char* s = expr_to_string(e);
    printf("%-14s %s\n", label, s);
    free(s);
}

int main(void) {
    // f = (x^2 - 1)/(x - 1) + 1/x, built from constructors (each takes its arguments)
    Expr* x2 = create_pow(create_var("x"), create_const(2));
    Expr* f = create_add(create_div(create_sub(x2, create_const(1)), create_sub(create_var("x"), create_const(1))),
                         create_div(create_const(1), create_var("x")));
    show("f =", f);
    Expr* s = simplify(f);                         // a reduced rational function
    show("simplified:", s);
    Expr* d = differentiate(f, "x");
    show("f' =", d);
    bool ok = true;
    printf("%-14s %.12g\n", "f(2.5) =", expr_eval(f, 2.5, &ok));

    // with exp: canonical Laurent form in exp(x) through simplify_exp
    const char* why = NULL;
    Expr* g = parse_input("(exp(2x) - 1)/(exp(x) - 1)");
    Expr* gs = simplify_exp(g, &why);
    show("simplify_exp:", gs);
    free_expr(f); free_expr(s); free_expr(d); free_expr(g); free_expr(gs);
    return 0;
}
```

Output:

```text
f =            (x^2 - 1) / (x - 1) + 1 / x
simplified:    (x^2 + x + 1) / (x)
f' =           (x^2 - 1) / (x^2)
f(2.5) =       3.9
simplify_exp:  1 + exp(x)
```

### Indefinite integrals

- **`integrate(f, "x")`** handles rational functions. It returns `NULL` on failure, with the reason in `integrate_error()`.
- **`integrate_exp(f)`** handles anything with exp, ln, sin, cos or tan. Its `RischResult` has one of three statuses:
  - `RISCH_OK` with the antiderivative;
  - `RISCH_NONELEMENTARY` with the part that has no elementary integral (a proof);
  - `RISCH_UNSUPPORTED` with a reason.
- **`integrate_rootsum`** works like `integrate` but always keeps a `RootSum` instead of rewriting it in radicals.

```c
// examples/integrate.c
// Indefinite integrals: rational functions (integrate), and exp/ln/sin/cos (integrate_exp,
// the Risch algorithm), which also proves when no elementary antiderivative exists.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static void rational(const char* text) {
    Expr* f = parse_input(text);
    Expr* F = integrate(f, "x");                    // NULL with integrate_error() on failure
    char* s = F ? expr_to_string(F) : NULL;
    printf("∫ %-22s = %s\n", text, s ? s : integrate_error());
    free(s); free_expr(F); free_expr(f);
}

static void risch(const char* text) {
    Expr* f = parse_input(text);
    RischResult r = integrate_exp(f);
    char* s = r.result ? expr_to_string(r.result) : NULL;
    if (r.status == RISCH_OK) printf("∫ %-22s = %s\n", text, s);
    else if (r.status == RISCH_NONELEMENTARY) printf("∫ %-22s : no elementary antiderivative (the part %s)\n", text, s);
    else printf("∫ %-22s : unsupported (%s)\n", text, r.reason);
    free(s); free_expr(r.result); free_expr(f);
}

int main(void) {
    rational("1/((x-1)(x^2+1))");
    rational("x/(x^4+1)");
    rational("1/(x^3-2)");                        // radicals: the roots of x^3 - 2 are known
    rational("1/(x^5+x+3)");                      // no radicals exist: an exact RootSum
    risch("x^2*exp(3x)");
    risch("x^2*ln(x)^3");
    risch("exp(x)/(1+exp(2x))");
    risch("x*exp(x)*sin(x)");
    risch("1/(2+cos(x))");
    risch("exp(x)/(exp(3x)+2)");                  // cubic residues: a RootSum over exp
    risch("exp(x^2)");                            // proved non-elementary
    risch("ln(x)/(x+1)");

    const char* why = NULL;                       // differentiate_exp: the same expressions
    Expr* g = parse_input("x*exp(x)*sin(x)");
    Expr* dg = differentiate_exp(g, &why);
    char* s = expr_to_string(dg);
    printf("d/dx x*exp(x)*sin(x) = %s\n", s);
    free(s); free_expr(dg); free_expr(g);
    return 0;
}
```

Output:

```text
∫ 1/((x-1)(x^2+1))       = 1/2 * ln(abs(x - 1)) - 1/4 * ln(abs(x^2 + 1)) - 1/2 * atan(x)
∫ x/(x^4+1)              = 1/2 * atan(x^2)
∫ 1/(x^3-2)              = 1/6 * 2^(1/3) * ln(abs(x - 2^(1/3))) - 1/12 * 2^(1/3) * ln(abs(x^2 + 2^(1/3) * x + 2^(2/3))) - 1/6 * 2^(1/3) * sqrt(3) * atan(1/3 * 2^(2/3) * sqrt(3) * x + 1/3 * sqrt(3))
∫ 1/(x^5+x+3)            = RootSum(253381t^5 - 160t^3 - 80t^2 - 15t - 1, t -> t * ln(x + (-16216384/16875t^4 + 4054096/16875t^3 - 111476/1875t^2 + 255941/16875t + 256/16875)))
∫ x^2*exp(3x)            = (1/3x^2 - 2/9x + 2/27) * exp(3x)
∫ x^2*ln(x)^3            = 1/3x^3 * ln(abs(x))^3 - 1/3x^3 * ln(abs(x))^2 + 2/9x^3 * ln(abs(x)) - 2/27x^3
∫ exp(x)/(1+exp(2x))     = atan(exp(x))
∫ x*exp(x)*sin(x)        = 1/2x * exp(x) * sin(x) + (-1/2x + 1/2) * exp(x) * cos(x)
∫ 1/(2+cos(x))           = 2/3 * sqrt(3) * atan(1/3 * sqrt(3) * tan(1/2x))
∫ exp(x)/(exp(3x)+2)     = RootSum(108t^3 - 1, t -> t * ln(exp(x) + 6 * t))
∫ exp(x^2)               : no elementary antiderivative (the part exp(x^2))
∫ ln(x)/(x+1)            : no elementary antiderivative (the part ln(abs(x)) / (x + 1))
d/dx x*exp(x)*sin(x) = (x + 1) * exp(x) * sin(x) + x * exp(x) * cos(x)
```

### Definite integrals

`integrate_definite` returns the exact value (or `NULL` with a `note` when only the number exists), the quadrature value `approx`, or an `error` (a pole, divergence). Bounds come from `parse_endpoint`, which reads `1/2`, `pi/4`, `inf` and `-inf`.

```c
// examples/definite.c
// Definite integrals: exact values (constants folded) checked against quadrature at run time.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static void definite(const char* text, const char* a_text, const char* b_text) {
    Rational a, b;
    rat_init(&a); rat_init(&b);
    Bounds ab = { { &a, 0, false }, { &b, 0, false } };
    Expr* f = parse_integrand(text);
    parse_endpoint(&a, &ab.a, a_text);            // "1/2", "pi/4", "inf", "-inf"
    parse_endpoint(&b, &ab.b, b_text);
    DefiniteResult r = integrate_definite(f, &ab);
    printf("∫ %-15s on [%s, %s]: ", text, a_text, b_text);
    if (r.error) printf("error: %s\n", r.error);
    else if (r.exact) { char* s = expr_to_string(r.exact); printf("= %s ~ %.12g\n", s, r.approx); free(s); }
    else printf("~ %.12g (%s)\n", r.approx, r.note);
    definite_free(&r);
    free_expr(f);
    rat_clear(&a); rat_clear(&b);
}

int main(void) {
    definite("1/(x^2+x+1)", "0", "1");
    definite("1/(x^4+1)", "-inf", "inf");
    definite("x^2*sin(x)", "0", "pi");
    definite("exp(-x)*cos(x)", "0", "inf");
    definite("1/(2+cos(x))", "0", "2pi");          // the jump of atan(tan(x/2)/√3) is corrected
    definite("exp(-x^2)", "-inf", "inf");          // no elementary antiderivative: the number
    definite("1/x", "-1", "1");                    // a pole in range
    return 0;
}
```

Output:

```text
∫ 1/(x^2+x+1)     on [0, 1]: = 1/9 * sqrt(3) * pi ~ 0.604599788078
∫ 1/(x^4+1)       on [-inf, inf]: = 1/2 * sqrt(2) * pi ~ 2.22144146908
∫ x^2*sin(x)      on [0, pi]: = pi^2 - 4 ~ 5.86960440109
∫ exp(-x)*cos(x)  on [0, inf]: = 1/2 ~ 0.499999999999
∫ 1/(2+cos(x))    on [0, 2pi]: = 2/3 * sqrt(3) * pi ~ 3.62759872847
∫ exp(-x^2)       on [-inf, inf]: ~ 1.77245385091 (no elementary antiderivative, so only the numeric value is given)
∫ 1/x             on [-1, 1]: error: the integrand has a pole in [a, b], so the integral diverges
```

### Linear algebra

`MAT_AT(m, r, c)` is the `Rational` at row r, column c. `linsys_solve` reports `LINSYS_UNIQUE`, `LINSYS_INCONSISTENT` or `LINSYS_INFINITE` along with the rank.

```c
// examples/linear.c
// Exact linear algebra: matrices of rationals, parsed from bracket notation or filled directly.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, char* s) {
    printf("%-10s %s\n", label, s);
    free(s);
}

int main(void) {
    Matrix* H = matrix_create(3, 3);              // the Hilbert matrix 1/(i + j + 1)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) rat_set_si(&MAT_AT(H, i, j), 1, i + j + 1);
    show("H =", matrix_to_string(H));
    Rational d;
    rat_init(&d);
    matrix_det(&d, H);
    show("det H =", rat_to_string(&d));
    Matrix* Hi = matrix_inverse(H);
    show("H^-1 =", matrix_to_string(Hi));

    Matrix* A = matrix_parse("[[1, 2, 3], [4, 5, 6], [7, 8, 10]]");
    Vector* b = vector_parse("[1, 0.5, 1/3]");
    LinSysResult* x = linsys_solve(A, b);
    show("A x = b:", linsys_result_to_string(x));
    Matrix* S = matrix_parse("[[1, 2, 3], [2, 4, 6], [1, 0, 1]]");
    printf("%-10s %d\n", "rank S =", matrix_rank(S));
    Matrix* R = matrix_rref(S);
    show("rref S =", matrix_to_string(R));

    rat_clear(&d);
    matrix_free(H); matrix_free(Hi); matrix_free(A); matrix_free(S); matrix_free(R);
    vector_free(b); linsys_result_free(x);
    return 0;
}
```

Output:

```text
H =        [[1, 1/2, 1/3]
 [1/2, 1/3, 1/4]
 [1/3, 1/4, 1/5]]
det H =    1/2160
H^-1 =     [[9, -36, 30]
 [-36, 192, -180]
 [30, -180, 180]]
A x = b:   x1 = -1, x2 = 1/2, x3 = 1/3
rank S =   2
rref S =   [[1, 0, 1]
 [0, 1, 1]
 [0, 0, 0]]
```

### Big numbers

```c
// examples/big_numbers.c
// Big numbers: every coefficient is an arbitrary-precision GMP rational, so nothing overflows
// or rounds, in factoring, roots, integrals or matrices.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static Expr* simplified(const char* text) {
    Expr* raw = parse_input(text);
    Expr* e = canonicalize(raw);
    free_expr(raw);
    return e;
}

int main(void) {
    // (x - 10^30)(x + 2^100): both roots come back exactly
    Expr* p = simplified("(x - 1000000000000000000000000000000)(x + 1267650600228229401496703205376)");
    char* ps = expr_to_string(p);
    printf("p        = %s\n", ps);
    Factorization* f = factorize(p->poly);
    char* fs = factorization_to_string(f);
    printf("factor p = %s\n", fs);
    Solutions* sol = solve_polynomial(p->poly);
    char* ss = solutions_to_string(sol);
    printf("solve p  = %s\n", ss);

    Expr* q = simplified("(x + 1)^60");           // binomial coefficients up to C(60, 30)
    char* qs = expr_to_string(q);
    printf("(x+1)^60 = %.120s...\n", qs);

    Expr* r = parse_input("1/(x^2 + 123456789012345678901234567890)");
    Expr* R = integrate(r, "x");
    char* Rs = expr_to_string(R);
    printf("∫ 1/(x^2 + 123456789012345678901234567890) = %s\n", Rs);

    Matrix* m = matrix_parse("[[1e40, 3], [7, 1/99999999999999999999]]");
    Rational d;
    rat_init(&d);
    matrix_det(&d, m);
    char* ds = rat_to_string(&d);
    printf("det      = %s\n", ds);

    free(ps); free(fs); free(ss); free(qs); free(Rs); free(ds);
    rat_clear(&d);
    matrix_free(m);
    free_factorization(f); free_solutions(sol);
    free_expr(p); free_expr(q); free_expr(r); free_expr(R);
    return 0;
}
```

Output:

```text
p        = x^2 + 267650600228229401496703205376x - 1267650600228229401496703205376000000000000000000000000000000
factor p = (x + 1267650600228229401496703205376)(x - 1000000000000000000000000000000)
solve p  = x = -1267650600228229401496703205376, x = 1000000000000000000000000000000
(x+1)^60 = x^60 + 60x^59 + 1770x^58 + 34220x^57 + 487635x^56 + 5461512x^55 + 50063860x^54 + 386206920x^53 + 2558620845x^52 + 147831...
∫ 1/(x^2 + 123456789012345678901234567890) = 1/41152263004115226300411522630 * sqrt(13717421001371742100137174210) * atan(1/3x / sqrt(13717421001371742100137174210))
det      = 9999999999999999997900000000000000000021/99999999999999999999
```

## Expression syntax

- Operators `+`, `-`, `*`, `/`, `^`, with parentheses and implicit multiplication (`5x^3`, `(x + 1)(x - 1)`, `2exp(x)`).
- Unary minus: `-x^2` means `-(x^2)`.
- Number literals are integers, decimals and exponents (`12`, `0.75`, `2.5e-3`), all exact.
- Exponents must be integers: `x^(1/2)` is rejected, and so is division by zero.
- Anything left over after a complete expression (e.g. `= 4`) is an error.
- Functions:
  - `exp(u)`, `ln(u)`, `sin(u)`, `cos(u)`, and `tan(u)` (read as `sin(u)/cos(u)`, printed back as `tan(u)`);
  - accepted by `simplify`, `differentiate` and `integrate`, with the limits below;
  - `factor` and `solve` take polynomials only.

## Limits

- When a logarithm's coefficient is an algebraic number of degree 3 or more, the answer is an exact `RootSum(q(t), t -> t * ln(S(t, x)))`: the sum of `t·ln(S(t, x))` over the roots of `q`. For rational functions this happens only when the roots aren't known in closed form. Over exp, ln or tan (and an outer ln) the RootSum is always kept, with S a polynomial in the tower written as `Σ t^i·(…)`; definite integrals over it give the numeric value.
- Complex roots of other irreducible factors of degree 3 or more (for example `x^3 + x + 1`, or the Capelli factors `x^4 ± 2x^2 + 2` of `x^8 + 4`) are counted but not located.
- Exp, sin and cos arguments must be polynomials. A denominator involving exp may use only one exponential tower: every `exp` argument a rational multiple of one polynomial. A denominator involving sin or cos needs every argument to be an integer multiple of one polynomial (`sin(x)` and `cos(2x)`). These are reported as not supported yet:
  - denominators whose residues are algebraic numbers of degree 7 or more, such as `exp(x)/(exp(7x)+2)` or `1/(sin(x)^7+2)`. These are elementary;
  - `exp` or `ln` together with a sin/cos denominator, such as `exp(x)/(1+cos(x))`;
  - `simplify` of any sum of sin/cos in a denominator (`tan(x)`), and `differentiate` of those that need `tan(x/2)`, such as `1/(1+cos(x))` (the exact derivative would be unreadable in `tan(x/2)`);
  - denominators with independent exponentials, such as `exp(x)/(exp(x)+exp(x^2))`;
  - `exp(1/x)`, `sin(1/x)`, and constants such as `exp(3)` or `cos(1)`;
  - `ln(2x)·ln(x)` and `exp(x)·ln(exp(x))`, where the outer ln lies in the inner tower (a new constant `ln 2`, or `ln(exp(x)) = x`);
  - an outer ln in a denominator whose residues are algebraic numbers of degree 7 or more, or of degree above 6 in the outer ln, three or more levels, `exp` of a non-polynomial (`exp(exp(x))`), and `ln` together with `sin` or `cos`;
  - definite integrals with `ln`; with a sum of sin/cos or exp in a denominator over an infinite range, or whose denominator has coefficients depending on x (`1/(x+exp(x))`) or a nonlinear argument (`1/(2+cos(x^2))`); and where such a denominator vanishes at a point where the numerator may too (`x/sin(x)` on `[-1, 1]`, which is finite there);
  - definite integrals over a range where a coefficient's denominator vanishes, even when the integrand stays finite (`sin(x)/x` on `[0, 1]`);
  - infinite ranges where some exp/sin/cos term does not decay exponentially: `cos(x)/x` on `[1, inf]` converges but is refused, and so is the divergent `exp(x)` on `[0, inf]`.
- Definite integrals: when the antiderivative has a `RootSum`, only the numeric value is given. Exact values are simplified in one quadratic field Q(√m): logs merge, atan pairs combine (`∫₀¹ 1/(x⁴+1)` gives `1/4 * sqrt(2) * ln(1 + sqrt(2)) + 1/8 * sqrt(2) * pi`). Values that need two radicals (`cos(pi/12)`) or higher roots (`2^(1/4)`) are folded but not merged.
- Decimals beside exact answers are doubles (10 significant digits); the exact part is the answer. The run-time quadrature check also works in doubles, so a definite integral whose value is beyond about 1e308 (`x^2` on `[0, 1e200]`) is given exactly with a note that it was not compared with quadrature.

## Tests

```bash
make && ./test_symbolic_solver
```

- 96 unit tests covering rational arithmetic, polynomial algebra, parsing, canonicalization, factoring, solving, calculus and linear systems (including Hilbert matrices).
- Regression tests for each previously found bug: big literals, truncated output, dropped integral terms, lost factor content, `-x^2` parsed as `(-x)^2`, trailing input silently ignored, high-degree and huge-range definite integrals, and roots near zero.
- Seeded property tests that check invariants on generated inputs:
  - `content · Π factors == input` exactly
  - each root interval brackets a sign change
  - the exact derivative of every antiderivative equals the integrand (for `RootSum` terms, compared exactly at rational points via traces over Q(α)); answers with radicals and trig constants are checked numerically in double precision, and to 60 digits by the SymPy cross-check
  - every exact root printed for a binomial/cyclotomic factor satisfies the factor
  - every exact definite integral agrees with adaptive quadrature (the run-time check never fails), including infinite bounds on generated convergent integrands
  - `A·solve(A, A·x) == A·x`
  - parse(print(p)) == p

  Set `POLY_SEED=<n>` to try other seeds.

`make check` is the quick lane after a change, a few seconds in all. It runs the suite, which prints its CPU time, and then replays `golden.jsonl`: about 1,900 CLI calls taken from the SymPy cross-check's generators, each of which must give exactly the recorded output. `POLY_TIME=1` prints the CPU time of each regression and property test.

A differential check against SymPy is also available. It is optional and needs `pip install sympy`:

```bash
make cli && python3 crosscheck_sympy.py --trials 40   # a few minutes; --seed N for others
```

It confirms each factorization multiplies back, each root and root count (every exactly printed root is substituted back at 50 digits), each antiderivative by differentiating it, and each definite integral against `mpmath.quad` at 40 digits. A separate stream uses infinite bounds and also checks that divergence is reported exactly when there is a real pole in range or the degree gap is below 2. Answers containing `RootSum`, radicals or trig constants are checked at 60 significant digits, including on a separate stream of binomial and cyclotomic denominators. The Risch streams compare every verdict with SymPy's `risch_integrate`.

## Architecture

1. **Rational arithmetic** (`rational_num.c`): GMP `mpq_t`, exact literal parsing, decimal output
2. **Polynomial algebra** (`polynomial.c`): sparse sorted terms built in order, dense multiply and divide, Euclidean GCD/XGCD, exact Horner evaluation
3. **Modular and integer polynomials** (`poly_zp.c`, `poly_zx.c`): dense F_p[x] and Z[x]
4. **Expression tree** (`expression.c`): tagged union with canonicalization to polynomial/rational form
5. **Parser** (`token.c`, `parser.c`): shunting-yard with implicit multiplication and unary minus
6. **Factoring** (`factor.c`, `zx_factor.c`, `modgcd.c`): Yun square-free decomposition (modular gcds, exactly checked) + Berlekamp–Zassenhaus
7. **Roots and solver** (`roots.c`, `solver.c`): Descartes/VCA isolation, exact refinement, radicals
8. **Integration** (`integrate.c`, `logpart.c`, `subres.c`, `rootform.c`, `surd.c`, `partial_frac.c`): Hermite reduction, log/atan terms, Lazard–Rioboo–Trager, radical forms, partial fractions
9. **Risch** (`expfield.c`, `risch.c`, `kpoly.c`, `residue.c`, `logpair_tmpl.h`, `rootsum_tmpl.h`, `cfield.c`, `expsum.c`, `tanfield.c`, `risch2.c`, `k2poly.c`): exp/ln towers, residues, exponential sums, tan(x) and tan(x/2), two levels
10. **Definite integrals** (`definite.c`, `deftower.c`, `constfold.c`): pole detection, exact evaluation with atan jump corrections, constant folding with Q(√m) arithmetic (`qsurd.c`), quadrature check
11. **Linear algebra** (`matrix.c`): Gaussian elimination over exact rationals

See [INTERNALS.md](INTERNALS.md) for how each algorithm works.

## License

MIT
