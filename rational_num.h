#ifndef RATIONAL_NUM_H
#define RATIONAL_NUM_H

#include <gmp.h>
#include <stdbool.h>

typedef struct {
    mpq_t val;
} Rational;

// lifecycle — every Rational must be initialized before use and cleared after
void rat_init(Rational* r);
void rat_clear(Rational* r);
void rat_init_set(Rational* dst, const Rational* src);

// construction (initializes and sets)
void rat_set_si(Rational* r, long num, long den);
void rat_set(Rational* dst, const Rational* src);

// convenience constructors that init + set in one call
Rational rat_from_si(long num, long den);
Rational rat_from_int(long n);
Rational rat_zero_val(void);
Rational rat_one_val(void);

// arithmetic — result must be pre-initialized
void rat_add(Rational* result, const Rational* a, const Rational* b);
void rat_sub(Rational* result, const Rational* a, const Rational* b);
void rat_mul(Rational* result, const Rational* a, const Rational* b);
void rat_div(Rational* result, const Rational* a, const Rational* b);
void rat_neg(Rational* result, const Rational* a);
void rat_abs(Rational* result, const Rational* a);
void rat_pow_int(Rational* result, const Rational* base, int exp);

// comparison
bool rat_eq(const Rational* a, const Rational* b);
bool rat_is_zero(const Rational* a);
bool rat_is_one(const Rational* a);
bool rat_is_negative(const Rational* a);
int  rat_cmp(const Rational* a, const Rational* b);

// conversion
double rat_to_double(const Rational* r);
char*  rat_to_string(const Rational* r); // caller frees

// integer GCD (for poly_content etc.)
void gcd_int(mpz_t result, const mpz_t a, const mpz_t b);

// access numerator/denominator as long (for small values)
long rat_num_si(const Rational* r);
long rat_den_si(const Rational* r);

#endif
