#ifndef RATIONAL_NUM_H
#define RATIONAL_NUM_H

#include <stdbool.h>

// typedef for future GMP swap
typedef long long rat_int_t;

typedef struct {
    rat_int_t num;
    rat_int_t den; // always > 0 after normalization
} Rational;

// construction
Rational rat_create(rat_int_t num, rat_int_t den);
Rational rat_from_int(rat_int_t n);
Rational rat_zero(void);
Rational rat_one(void);

// arithmetic — all return new values, inputs unmodified
Rational rat_add(Rational a, Rational b);
Rational rat_sub(Rational a, Rational b);
Rational rat_mul(Rational a, Rational b);
Rational rat_div(Rational a, Rational b);
Rational rat_neg(Rational a);
Rational rat_abs(Rational a);
Rational rat_pow_int(Rational base, int exp);

// comparison
bool rat_eq(Rational a, Rational b);
bool rat_is_zero(Rational a);
bool rat_is_one(Rational a);
bool rat_is_negative(Rational a);
int  rat_cmp(Rational a, Rational b);

// conversion
double rat_to_double(Rational r);
char*  rat_to_string(Rational r); // caller frees

// internal
rat_int_t gcd_int(rat_int_t a, rat_int_t b);

#endif
