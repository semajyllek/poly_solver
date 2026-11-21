# Changelog

## 0.1.0

The first packaged release.

- **Library:** one public header, `polysolver.h`, with the whole API prefixed `ps_` (`PS_` for
  macros and enum constants). The static library defines only `ps_` and `psi_` symbols, and the
  shared library exports only `ps_` ones, so programs can use any other names.
- **Build:** a portable Makefile (GMP found with pkg-config, `PREFIX`/`DESTDIR`, `install`,
  `uninstall`, `dist`) and CMake (`find_package(polysolver)`, `polysolver::polysolver`), both
  installing the header, static and shared libraries, the CLI and `polysolver.pc`.
- **Layout:** `src/`, `include/`, `cli/`, `tests/`, `tools/`, `examples/`.
- **Engine:**
  - simplification, factoring over Q, exact and certified real roots;
  - linear algebra over Q;
  - integration:
    - rational functions (Hermite, Lazard–Rioboo–Trager, closed-form RootSums);
    - exp, ln, sin, cos and tan through the Risch algorithm, with proofs of non-elementarity;
    - one square root of a quadratic, including conics without a rational point;
  - definite integrals with exact limits, atan jump corrections and a quadrature check, on
    rational, π-multiple and infinite bounds.
