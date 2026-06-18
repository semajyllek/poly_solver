#include "rational_num.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

// Euclidean GCD on absolute values
rat_int_t gcd_int(rat_int_t a, rat_int_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        rat_int_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

static Rational normalize(rat_int_t num, rat_int_t den) {
    assert(den != 0);
    if (num == 0) return (Rational){0, 1};
    // sign always in numerator
    if (den < 0) { num = -num; den = -den; }
    rat_int_t g = gcd_int(num, den);
    return (Rational){num / g, den / g};
}

Rational rat_create(rat_int_t num, rat_int_t den) {
    return normalize(num, den);
}

Rational rat_from_int(rat_int_t n) {
    return (Rational){n, 1};
}

Rational rat_zero(void) {
    return (Rational){0, 1};
}

Rational rat_one(void) {
    return (Rational){1, 1};
}

Rational rat_add(Rational a, Rational b) {
    // a.num/a.den + b.num/b.den = (a.num*b.den + b.num*a.den) / (a.den*b.den)
    // reduce intermediate overflow risk by dividing by gcd of denominators first
    rat_int_t g = gcd_int(a.den, b.den);
    rat_int_t lcm_den = a.den / g * b.den;
    rat_int_t num = a.num * (b.den / g) + b.num * (a.den / g);
    return normalize(num, lcm_den);
}

Rational rat_sub(Rational a, Rational b) {
    b.num = -b.num;
    return rat_add(a, b);
}

Rational rat_mul(Rational a, Rational b) {
    // cross-reduce before multiplying to limit overflow
    rat_int_t g1 = gcd_int(a.num, b.den);
    rat_int_t g2 = gcd_int(b.num, a.den);
    rat_int_t num = (a.num / g1) * (b.num / g2);
    rat_int_t den = (a.den / g2) * (b.den / g1);
    return normalize(num, den);
}

Rational rat_div(Rational a, Rational b) {
    assert(b.num != 0);
    return rat_mul(a, (Rational){b.den, b.num});
}

Rational rat_neg(Rational a) {
    return (Rational){-a.num, a.den};
}

Rational rat_abs(Rational a) {
    return (Rational){a.num < 0 ? -a.num : a.num, a.den};
}

Rational rat_pow_int(Rational base, int exp) {
    if (exp < 0) {
        assert(base.num != 0);
        base = (Rational){base.den, base.num};
        if (base.den < 0) { base.num = -base.num; base.den = -base.den; }
        exp = -exp;
    }
    Rational result = rat_one();
    while (exp > 0) {
        if (exp & 1) result = rat_mul(result, base);
        base = rat_mul(base, base);
        exp >>= 1;
    }
    return result;
}

bool rat_eq(Rational a, Rational b) {
    return a.num == b.num && a.den == b.den;
}

bool rat_is_zero(Rational a) {
    return a.num == 0;
}

bool rat_is_one(Rational a) {
    return a.num == 1 && a.den == 1;
}

bool rat_is_negative(Rational a) {
    return a.num < 0;
}

int rat_cmp(Rational a, Rational b) {
    // a/b vs c/d => compare a*d vs c*b (denominators always positive)
    rat_int_t lhs = a.num * b.den;
    rat_int_t rhs = b.num * a.den;
    if (lhs < rhs) return -1;
    if (lhs > rhs) return 1;
    return 0;
}

double rat_to_double(Rational r) {
    return (double)r.num / (double)r.den;
}

char* rat_to_string(Rational r) {
    char buf[64];
    if (r.den == 1)
        snprintf(buf, sizeof(buf), "%lld", r.num);
    else
        snprintf(buf, sizeof(buf), "%lld/%lld", r.num, r.den);

    char* s = malloc(strlen(buf) + 1);
    strcpy(s, buf);
    return s;
}
