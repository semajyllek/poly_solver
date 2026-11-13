#define _POSIX_C_SOURCE 200809L
#include "rational_num.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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

// --- exact literal parsing / decimal output ---

// read a run of digits starting at s into an mpz; returns pointer past the run
static const char* read_digits(mpz_t out, const char* s) {
    const char* e = s;
    while (isdigit((unsigned char)*e)) e++;
    char* tmp = strndup(s, (size_t)(e - s));
    mpz_set_str(out, tmp, 10);
    free(tmp);
    return e;
}

// multiply r by 10^k (k may be negative)
static void rat_scale_pow10(Rational* r, long k) {
    mpz_t p;
    mpz_init(p);
    mpz_ui_pow_ui(p, 10, (unsigned long)(k < 0 ? -k : k));
    if (k >= 0) mpz_mul(mpq_numref(r->val), mpq_numref(r->val), p);
    else mpz_mul(mpq_denref(r->val), mpq_denref(r->val), p);
    mpq_canonicalize(r->val);
    mpz_clear(p);
}

const char* rat_parse_literal(Rational* out, const char* s) {
    if (!isdigit((unsigned char)*s)) return NULL;
    mpz_t whole;
    mpz_init(whole);
    const char* p = read_digits(whole, s);
    long scale = 0; // value = whole * 10^scale
    if (*p == '.' && isdigit((unsigned char)p[1])) {
        const char* frac = p + 1;
        const char* fe = frac;
        while (isdigit((unsigned char)*fe)) fe++;
        // concatenate integer and fractional digits, then shift
        size_t n1 = (size_t)(p - s), n2 = (size_t)(fe - frac);
        char* all = malloc(n1 + n2 + 1);
        memcpy(all, s, n1);
        memcpy(all + n1, frac, n2);
        all[n1 + n2] = '\0';
        mpz_set_str(whole, all, 10);
        free(all);
        scale = -(long)n2;
        p = fe;
    }
    // optional exponent, only when digits follow (so "2e" stays 2 * e)
    if ((*p == 'e' || *p == 'E')) {
        const char* q = p + 1;
        int sign = 1;
        if (*q == '+' || *q == '-') { sign = (*q == '-') ? -1 : 1; q++; }
        if (isdigit((unsigned char)*q)) {
            long ev = 0;
            while (isdigit((unsigned char)*q) && ev < 100000) ev = ev * 10 + (*q++ - '0');
            while (isdigit((unsigned char)*q)) q++;
            scale += sign * ev;
            p = q;
        }
    }
    mpq_set_z(out->val, whole);
    rat_scale_pow10(out, scale);
    mpz_clear(whole);
    return p;
}

char* rat_to_decimal(const Rational* r, int digits) {
    mpz_t scaled, p10;
    mpz_init(scaled);
    mpz_init(p10);
    mpz_ui_pow_ui(p10, 10, (unsigned long)digits);
    mpz_mul(scaled, mpq_numref(r->val), p10);
    mpz_tdiv_q(scaled, scaled, mpq_denref(r->val)); // truncate toward zero
    bool neg = mpz_sgn(scaled) < 0 || (mpz_sgn(scaled) == 0 && mpq_sgn(r->val) < 0);
    mpz_abs(scaled, scaled);
    char* ds = mpz_get_str(NULL, 10, scaled);
    size_t len = strlen(ds);
    // left-pad so there is at least one digit before the point
    size_t width = len > (size_t)digits ? len : (size_t)digits + 1;
    char* out = malloc(width + 3);
    char* w = out;
    if (neg) *w++ = '-';
    size_t pad = width - len;
    for (size_t i = 0; i < width; i++) {
        *w++ = (i < pad) ? '0' : ds[i - pad];
        if (i == width - (size_t)digits - 1 && digits > 0) *w++ = '.';
    }
    *w = '\0';
    free(ds);
    mpz_clear(scaled);
    mpz_clear(p10);
    return out;
}

#define SQRT_TRIAL_LIMIT 10000UL

void rat_sqrt_split(const Rational* q, Rational* r, mpz_t m) {
    // sqrt(n/d) = sqrt(n·d) / d, then pull square factors out of n·d
    mpz_t s, p2;
    mpz_init_set_ui(s, 1);
    mpz_init(p2);
    mpz_mul(m, mpq_numref(q->val), mpq_denref(q->val));
    for (unsigned long p = 2; p < SQRT_TRIAL_LIMIT && mpz_cmp_ui(m, 1) > 0; p++) {
        mpz_set_ui(p2, p * p);
        while (mpz_divisible_p(m, p2)) {
            mpz_divexact(m, m, p2);
            mpz_mul_ui(s, s, p);
        }
    }
    if (mpz_perfect_square_p(m)) { // a large square cofactor that trial division missed
        mpz_sqrt(p2, m);
        mpz_mul(s, s, p2);
        mpz_set_ui(m, 1);
    }
    mpq_set_num(r->val, s);
    mpq_set_den(r->val, mpq_denref(q->val));
    mpq_canonicalize(r->val);
    mpz_clear(s);
    mpz_clear(p2);
}
