# Roadmap

Planned work, grouped by area. Within each group, items are roughly in the order they would be
done. Sizes are relative scope, not calendar time: **S** is a contained change, **M** is a new
component with its own tests, **L** spans several components or needs a new algorithmic layer.
When an item ships, it moves to [CHANGELOG.md](CHANGELOG.md).

## First: reach and trust

These do the most for people other than us: an honest public comparison, a stable release to
depend on, the languages most users work in, and integrals they can actually state.

1. **Integration benchmark (S).** Run a standard corpus, such as Rubi's published test suite of
   tens of thousands of integrals with known answers. Record the solve rate, wrong answers
   (there should be none) and time per family, beside SymPy (and FriCAS where available). This
   gives the "how good is it" claim a number and reorders the rest of this list.
2. **Confirm the public names, then tag v0.1.0 (S)** to publish the first release.
3. **Python bindings (M),** a thin package on PyPI. This is where most potential users are.
4. **WebAssembly build (M):** compile the library to WASM (Emscripten) with a small JavaScript
   API and an npm package, so web pages, notebooks and calculators get exact integration on the
   client. Nothing comparable exists for JavaScript.
5. **Symbolic parameters and several variables (L).** `∫ 1/(x^2+a^2)`, with case splits where
    the answer depends on a parameter's sign or on a parameter being zero. For users, this is
    the largest single gap.

## Measurement

6. **Polynomial benchmark against FLINT (S).** Multiplication, gcd, factoring (including
   Swinnerton-Dyer and other hard recombination cases) and real-root isolation, at degrees
   10–1000 and coefficient sizes up to thousands of bits, to know the real gap at each size.
7. **Track both benchmarks per release (S),** as numbers to compare, not pass/fail (run locally while CI
   minutes are scarce).

## Square roots

8. **Non-elementarity proofs by Chebyshev's theorem (S–M).** ∫ xᵐ·(a + b·xⁿ)ᵖ is elementary
   exactly when p, (m+1)/n or (m+1)/n + p is an integer. `1/sqrt(x^3+1)` gets a proof instead of
   "not supported", and the elementary cases (`x^2/sqrt(x^3+1)`, `x^3*sqrt(x^4+1)`) integrate
   through their rational substitution.
9. **Square roots of rational functions (M),** such as `sqrt(1/x)` and `sqrt((x+1)/(x-1))`.
10. **Two different square roots (L),** such as `sqrt(x)*sqrt(x+1)` and `1/(sqrt(x)+sqrt(x+1))`.
11. **Square roots combined with exp, ln, sin or cos (L).**
12. **Definite integrals with square roots at multiples of π (S).**
13. **Poles whose residues generate only part of the roots' field (M),** such as
   `1/((x^4+1)*sqrt(3-x^2))`.
14. **`sqrt(x^2)` as |x| (M),** handled piecewise.
15. **General algebraic functions (L):** Trager–Bronstein integration of any radical or
    algebraic function, the step that would match FriCAS's breadth.

## RootSums and exact constants

16. **Real forms for RootSums over t⁴ + p·t² + q (M).** These roots have closed radical forms.
    Most of the remaining number-only definite integrals would become exact, on the rational,
    square-root and exp/ln paths alike.
17. **Exact definite values for any RootSum antiderivative (M),** beyond those closed forms.
18. **Nested 0/0 atan arguments at a rational point (S),** e.g. `1/((x^2+2)*sqrt(4-3x^2))` on
    [0, 1/2], which is number-only today.
19. **Merging constants that need two radicals or higher roots (M),** such as `cos(pi/12)` and
    `2^(1/4)`.

## Exp, ln, sin and cos (Risch)

20. **Residues of degree 7 or more (M),** such as `exp(x)/(exp(7x)+2)` and `1/(sin(x)^7+2)`.
21. **exp or ln together with a sin/cos denominator (M),** such as `exp(x)/(1+cos(x))`; also ln
    with sin or cos.
22. **Independent exponentials in a denominator (L),** such as `exp(x)/(exp(x)+exp(x^2))`.
23. **Non-polynomial arguments (L),** such as `exp(1/x)`, `sin(1/x)` and `exp(exp(x))`; also
    constants like `exp(3)` and `cos(1)`.
24. **An outer ln that lies in the inner tower (M),** such as `ln(2x)·ln(x)` and
    `exp(x)·ln(exp(x))`.
25. **Three or more tower levels (L),** and degree above 6 in the outer ln.
26. **`simplify` of sin/cos sums in a denominator (M),** such as `tan(x)`, and readable
    derivatives of the `tan(x/2)` cases, such as `1/(1+cos(x))`.

## Expressiveness

27. **Special functions in answers (L):** erf, Ei, the logarithmic integral, polylogarithms, and
    elliptic integrals. Where the library now proves "non-elementary", it would also give the
    answer in these, with their derivatives, numeric evaluation and definite-integral values.

## Definite integrals

28. **Definite integrals with ln (M).**
29. **Integrands that stay finite where a coefficient's denominator vanishes (M),** such as
    `sin(x)/x` on [0, 1] and `x/sin(x)` on [−1, 1].
30. **Slowly converging infinite ranges (M),** such as `cos(x)/x` on [1, ∞), refused today.
31. **More exp and sin/cos denominators (M):** over infinite ranges, with x-dependent
    coefficients (`1/(x+exp(x))`), or with nonlinear arguments (`1/(2+cos(x^2))`).
32. **Quadrature checks for values beyond about 1e308 (S),** by scaling or extended
    precision, instead of a note.

## Polynomial core speed

Each of these is a known algorithm of bounded size; together they close most of the gap to
FLINT and PARI at large sizes. Integration inputs are small, so these matter mainly for use as a
polynomial library.

33. **Integer coefficients with one common denominator (M).** Store a polynomial in Q[x] as a
    dense array of integers with one shared denominator, instead of a linked list of rationals.
    That removes an allocation and a fraction reduction from every coefficient operation. This is
    a large constant factor everywhere, but it touches the core.
34. **Fast multiplication (M):** Karatsuba, then Kronecker substitution (pack the polynomial into
    one big integer and use GMP's FFT multiplication) for large sizes, with measured crossovers.
35. **Heuristic gcd and half-gcd (M):** evaluate at a large integer and take one integer gcd
    (most inputs finish there), with half-gcd for large degrees and 64-bit primes in the modular
    gcd.
36. **Quadratic multifactor Hensel lifting (M):** doubling the precision each step, all factors
    at once through a factor tree, instead of one power of p at a time.
37. **Van Hoeij recombination (L):** lattice reduction (LLL) instead of subset search, which
    removes the exponential worst case (polynomials irreducible over Q with many modular
    factors). It needs a solid LLL implementation.
38. **Faster real-root isolation (M–L):** asymptotically fast Taylor shifts, then a method that
    works mostly in interval arithmetic and raises the precision only where needed (ANewDsc
    style), keeping the exact certificate at the end.
39. **Locating complex roots (M)** of irreducible factors of degree 3 or more (`x^3+x+1`, the
    factors `x^4 ± 2x^2 + 2` of `x^8+4`), with certified complex isolating boxes.
40. **Fraction-free and multimodular linear algebra (M):** Bareiss elimination, determinants by
    CRT, and Dixon's p-adic lifting for linear systems.
41. **Optional FLINT backend (M).** When FLINT is installed, hand large factoring, gcd and root
    isolation to it. FLINT is LGPL, which is compatible with this MIT library when FLINT is
    linked dynamically.

## Packaging and reach

42. **The first CI run (S),** which is also the first Linux build, when there are Actions minutes
    for it. The golden decimals may need a tolerance for Linux's libm.
43. **Homebrew tap (S)** for `brew install polysolver`.
44. **Debian/Ubuntu package (S–M).**
45. **API reference (S):** generated from the header's comments.

## Code health

46. **Delete the unused `src/token_list.c` and `src/token_stack.c` (S).**
47. **Share the atan-jump helpers (S)** between `definite.c` and `sqrtdef.c`.
48. **Bring the seven older functions with more than three parameters into line (S).**
49. **Keep the verification loop fast (ongoing).** Suite CPU is 2.24 s at -O2; profile before
    features, prefer algorithmic speedups, and keep `make check` to seconds.
