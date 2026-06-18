#include "rational_num.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// --- lifecycle ---

void rat_init(Rational* r) {
    mpq_init(r->val);
}

void rat_clear(Rational* r) {
    mpq_clear(r->val);
}

void rat_init_set(Rational* dst, const Rational* src) {
    mpq_init(dst->val);
    mpq_set(dst->val, src->val);
}

// --- construction ---

void rat_set_si(Rational* r, long num, long den) {
    mpq_set_si(r->val, num, (unsigned long)(den < 0 ? -den : den));
    if (den < 0) mpq_neg(r->val, r->val);
    mpq_canonicalize(r->val);
}

void rat_set(Rational* dst, const Rational* src) {
    mpq_set(dst->val, src->val);
}

Rational rat_from_si(long num, long den) {
    Rational r;
    mpq_init(r.val);
    rat_set_si(&r, num, den);
    return r;
}

Rational rat_from_int(long n) {
    return rat_from_si(n, 1);
}

Rational rat_zero_val(void) {
    Rational r;
    mpq_init(r.val);
    return r;
}

Rational rat_one_val(void) {
    return rat_from_int(1);
}

// --- arithmetic ---

void rat_add(Rational* result, const Rational* a, const Rational* b) {
    mpq_add(result->val, a->val, b->val);
}

void rat_sub(Rational* result, const Rational* a, const Rational* b) {
    mpq_sub(result->val, a->val, b->val);
}

void rat_mul(Rational* result, const Rational* a, const Rational* b) {
    mpq_mul(result->val, a->val, b->val);
}

void rat_div(Rational* result, const Rational* a, const Rational* b) {
    mpq_div(result->val, a->val, b->val);
}

void rat_neg(Rational* result, const Rational* a) {
    mpq_neg(result->val, a->val);
}

void rat_abs(Rational* result, const Rational* a) {
    mpq_abs(result->val, a->val);
}

void rat_pow_int(Rational* result, const Rational* base, int exp) {
    if (exp < 0) {
        Rational inv;
        rat_init(&inv);
        mpq_inv(inv.val, base->val);
        rat_pow_int(result, &inv, -exp);
        rat_clear(&inv);
        return;
    }
    mpq_set_si(result->val, 1, 1);
    Rational b;
    rat_init_set(&b, base);
    while (exp > 0) {
        if (exp & 1) mpq_mul(result->val, result->val, b.val);
        mpq_mul(b.val, b.val, b.val);
        exp >>= 1;
    }
    rat_clear(&b);
}

// --- comparison ---

bool rat_eq(const Rational* a, const Rational* b) {
    return mpq_equal(a->val, b->val) != 0;
}

bool rat_is_zero(const Rational* a) {
    return mpq_sgn(a->val) == 0;
}

bool rat_is_one(const Rational* a) {
    return mpq_cmp_si(a->val, 1, 1) == 0;
}

bool rat_is_negative(const Rational* a) {
    return mpq_sgn(a->val) < 0;
}

int rat_cmp(const Rational* a, const Rational* b) {
    return mpq_cmp(a->val, b->val);
}

// --- conversion ---

double rat_to_double(const Rational* r) {
    return mpq_get_d(r->val);
}

char* rat_to_string(const Rational* r) {
    // check if denominator is 1
    if (mpz_cmp_si(mpq_denref(r->val), 1) == 0) {
        char* s = mpz_get_str(NULL, 10, mpq_numref(r->val));
        return s;
    }
    char* s = mpq_get_str(NULL, 10, r->val);
    return s;
}

// --- integer GCD ---

void gcd_int(mpz_t result, const mpz_t a, const mpz_t b) {
    mpz_gcd(result, a, b);
}

// --- accessors ---

long rat_num_si(const Rational* r) {
    return mpz_get_si(mpq_numref(r->val));
}

long rat_den_si(const Rational* r) {
    return mpz_get_si(mpq_denref(r->val));
}
