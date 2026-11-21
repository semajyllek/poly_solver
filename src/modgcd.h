#ifndef MODGCD_H
#define MODGCD_H

#include "polynomial.h"

// Exact results over Q computed modularly: images in F_p for many word-size primes,
// combined by CRT, rational reconstruction, then an exact check over Q. Each returns NULL
// when the check never passes within the prime budget; the caller then uses the
// (slower) Euclidean algorithm over Q. A non-NULL result is always exact.

// monic gcd(a, b); a and b nonzero
Polynomial* modular_gcd(const Polynomial* a, const Polynomial* b);

// s with s·a ≡ 1 (mod M), deg s < deg M; a must be invertible modulo M
Polynomial* modular_inverse(const Polynomial* a, const Polynomial* M);

// monic gcd(a, b) (either may be zero): modular when both have degree >= 6, where Euclid's
// coefficient growth dominates; Euclid (poly_gcd) when small or if the modular route gives up
Polynomial* poly_gcd_fast(const Polynomial* a, const Polynomial* b);

#endif
