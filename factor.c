#include "factor.h"
#include "modgcd.h"
#include "strbuf.h"
#include "poly_zx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// gcd(a/b, c/d) = gcd(a,c) / lcm(b,d)
Rational poly_content(const Polynomial* poly) {
    if (poly->head == NULL) return rat_zero_val();
    Rational content;
    rat_init(&content);
    rat_abs(&content, &poly->head->coeff);
    Term* t = poly->head->next;
    while (t) {
        Rational abs_c;
        rat_init(&abs_c);
        rat_abs(&abs_c, &t->coeff);
        // gcd of two rationals via GMP integer operations
        mpz_t num_gcd, den_lcm, g;
        mpz_init(num_gcd); mpz_init(den_lcm); mpz_init(g);
        mpz_gcd(num_gcd, mpq_numref(content.val), mpq_numref(abs_c.val));
        mpz_gcd(g, mpq_denref(content.val), mpq_denref(abs_c.val));
        mpz_divexact(den_lcm, mpq_denref(content.val), g);
        mpz_mul(den_lcm, den_lcm, mpq_denref(abs_c.val));
        mpq_set_num(content.val, num_gcd);
        mpq_set_den(content.val, den_lcm);
        mpq_canonicalize(content.val);
        mpz_clear(num_gcd); mpz_clear(den_lcm); mpz_clear(g);
        rat_clear(&abs_c);
        t = t->next;
    }
    return content;
}

Polynomial* poly_primitive_part(const Polynomial* poly) {
    Rational c = poly_content(poly);
    if (rat_is_zero(&c)) { rat_clear(&c); return poly_copy(poly); }
    Polynomial* result = create_polynomial();
    Term* t = poly->head;
    while (t) {
        Rational div_r;
        rat_init(&div_r);
        rat_div(&div_r, &t->coeff, &c);
        add_term(result, &div_r, t->exponent);
        rat_clear(&div_r);
        t = t->next;
    }
    rat_clear(&c);
    return result;
}

// --- square-free decomposition (Yun) ---

static void sqf_push(SquareFreeDecomp* d, Polynomial* part) {
    d->parts = realloc(d->parts, (size_t)(d->count + 1) * sizeof(Polynomial*));
    d->parts[d->count++] = part;
}

// b - c' (one step of Yun's recurrence)
static Polynomial* minus_derivative(const Polynomial* c, const Polynomial* b) {
    Polynomial* bd = poly_derivative(b);
    Polynomial* r = poly_subtract(c, bd);
    free_polynomial(bd);
    return r;
}

SquareFreeDecomp square_free_decompose(const Polynomial* f) {
    SquareFreeDecomp d = { NULL, 0 };
    if (polynomial_degree(f) <= 0) return d;
    Polynomial* fp = poly_derivative(f);
    Polynomial* a0 = poly_gcd_fast(f, fp);
    Polynomial* b = poly_divide(f, a0);
    Polynomial* c = poly_divide(fp, a0);
    Polynomial* dd = minus_derivative(c, b);
    free_polynomial(fp); free_polynomial(a0); free_polynomial(c);
    // each pass peels off the factors of multiplicity exactly i (Yun 1976)
    while (polynomial_degree(b) > 0) {
        Polynomial* ai = poly_gcd_fast(b, dd);
        Polynomial* nb = poly_divide(b, ai);
        Polynomial* nc = poly_divide(dd, ai);
        free_polynomial(b); free_polynomial(dd);
        b = nb;
        dd = minus_derivative(nc, b);
        free_polynomial(nc);
        sqf_push(&d, ai); // monic part of multiplicity d.count
    }
    free_polynomial(b); free_polynomial(dd);
    return d;
}

void free_square_free(SquareFreeDecomp* d) {
    for (int i = 0; i < d->count; i++) free_polynomial(d->parts[i]);
    free(d->parts);
    d->parts = NULL;
    d->count = 0;
}

// --- full factorization ---

static void fact_push(Factorization* f, Polynomial* p, int mult) {
    f->factors = realloc(f->factors, (size_t)(f->count + 1) * sizeof(Polynomial*));
    f->multiplicities = realloc(f->multiplicities, (size_t)(f->count + 1) * sizeof(int));
    f->factors[f->count] = p;
    f->multiplicities[f->count++] = mult;
}

// irreducible integer factors of one square-free part, appended with multiplicity
static void factor_part(Factorization* f, const Polynomial* part, int mult) {
    Polynomial* prim = poly_primitive_part(part); // integer coefficients, gcd 1
    ZxPoly z = zx_from_poly(prim);
    zx_make_primitive(&z);
    ZxFactorList L = zx_factor_squarefree(&z);
    for (int i = 0; i < L.count; i++) fact_push(f, zx_to_poly(&L.items[i]), mult);
    zx_factor_list_free(&L);
    zx_free(&z);
    free_polynomial(prim);
}

// order: degree, then linear factors by root, then coefficients from the top
static int factor_cmp(const Polynomial* a, const Polynomial* b) {
    int da = polynomial_degree(a), db = polynomial_degree(b);
    if (da != db) return da < db ? -1 : 1;
    Rational x, y;
    rat_init(&x); rat_init(&y);
    int cmp = 0;
    if (da == 1) { // compare roots -c0/c1
        Rational c0, c1;
        rat_init(&c0); rat_init(&c1);
        poly_coeff(&c0, a, 0); poly_coeff(&c1, a, 1); rat_div(&x, &c0, &c1); rat_neg(&x, &x);
        poly_coeff(&c0, b, 0); poly_coeff(&c1, b, 1); rat_div(&y, &c0, &c1); rat_neg(&y, &y);
        cmp = rat_cmp(&x, &y);
        rat_clear(&c0); rat_clear(&c1);
    }
    for (int e = da; e >= 0 && cmp == 0; e--) {
        poly_coeff(&x, a, e);
        poly_coeff(&y, b, e);
        cmp = rat_cmp(&x, &y);
    }
    rat_clear(&x); rat_clear(&y);
    return cmp;
}

// insertion sort of factors (with their multiplicities) for deterministic output
static void sort_factors(Factorization* f) {
    for (int i = 1; i < f->count; i++) {
        for (int j = i; j > 0 && factor_cmp(f->factors[j - 1], f->factors[j]) > 0; j--) {
            Polynomial* tp = f->factors[j]; f->factors[j] = f->factors[j - 1]; f->factors[j - 1] = tp;
            int tm = f->multiplicities[j]; f->multiplicities[j] = f->multiplicities[j - 1];
            f->multiplicities[j - 1] = tm;
        }
    }
}

// content = lc(poly) / lc(Π factors^m), so that poly == content · Π factors^m exactly
static void set_content(Factorization* f, const Polynomial* poly) {
    Rational lc_prod = rat_one_val(), lc_pow;
    rat_init(&lc_pow);
    for (int i = 0; i < f->count; i++) {
        rat_pow_int(&lc_pow, &f->factors[i]->head->coeff, f->multiplicities[i]);
        rat_mul(&lc_prod, &lc_prod, &lc_pow);
    }
    rat_div(&f->content, &poly->head->coeff, &lc_prod);
    rat_clear(&lc_prod);
    rat_clear(&lc_pow);
}

Factorization* factorize(const Polynomial* poly) {
    Factorization* f = malloc(sizeof(Factorization));
    f->factors = NULL;
    f->multiplicities = NULL;
    f->count = 0;
    f->content = rat_zero_val();
    if (poly->head == NULL) return f;
    SquareFreeDecomp d = square_free_decompose(poly);
    for (int i = 0; i < d.count; i++)
        if (polynomial_degree(d.parts[i]) > 0) factor_part(f, d.parts[i], i + 1);
    free_square_free(&d);
    sort_factors(f);
    set_content(f, poly);
    return f;
}

void free_factorization(Factorization* f) {
    if (!f) return;
    for (int i = 0; i < f->count; i++)
        free_polynomial(f->factors[i]);
    free(f->factors);
    free(f->multiplicities);
    rat_clear(&f->content);
    free(f);
}

char* factorization_to_string(const Factorization* f) {
    StrBuf sb = sb_new();
    Rational neg_one = rat_from_int(-1);
    bool unit = rat_is_one(&f->content) || rat_eq(&f->content, &neg_one);
    if (f->count == 0 || !unit) sb_append_owned(&sb, rat_to_string(&f->content));
    else if (rat_is_negative(&f->content)) sb_append(&sb, "-");
    rat_clear(&neg_one);
    for (int i = 0; i < f->count; i++) {
        sb_append(&sb, "(");
        sb_append_owned(&sb, poly_to_string(f->factors[i]));
        sb_append(&sb, ")");
        if (f->multiplicities[i] > 1) sb_appendf(&sb, "^%d", f->multiplicities[i]);
    }
    return sb_finish(&sb);
}
