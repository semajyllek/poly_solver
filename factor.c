#include "factor.h"
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

Polynomial* synthetic_divide(const Polynomial* poly, const Rational* root) {
    int deg = polynomial_degree(poly);
    if (deg <= 0) return create_polynomial();

    // dense coefficient array (descending order)
    Rational* coeffs = malloc((deg + 1) * sizeof(Rational));
    for (int i = 0; i <= deg; i++) coeffs[i] = rat_zero_val();
    Term* t = poly->head;
    while (t) {
        rat_set(&coeffs[deg - t->exponent], &t->coeff);
        t = t->next;
    }

    Rational* res = malloc(deg * sizeof(Rational));
    rat_init_set(&res[0], &coeffs[0]);
    for (int i = 1; i < deg; i++) {
        rat_init(&res[i]);
        Rational tmp;
        rat_init(&tmp);
        rat_mul(&tmp, &res[i - 1], root);
        rat_add(&res[i], &coeffs[i], &tmp);
        rat_clear(&tmp);
    }

    Polynomial* quotient = create_polynomial();
    for (int i = 0; i < deg; i++) {
        if (!rat_is_zero(&res[i]))
            add_term(quotient, &res[i], deg - 1 - i);
    }

    for (int i = 0; i <= deg; i++) rat_clear(&coeffs[i]);
    for (int i = 0; i < deg; i++) rat_clear(&res[i]);
    free(coeffs);
    free(res);
    return quotient;
}

// evaluate polynomial at a rational point (exact)
static void poly_eval_rational(Rational* result, const Polynomial* poly, const Rational* x) {
    if (poly->head == NULL) { mpq_set_si(result->val, 0, 1); return; }
    Term* t = poly->head;
    rat_set(result, &t->coeff);
    int prev_exp = t->exponent;
    t = t->next;
    while (t) {
        int gap = prev_exp - t->exponent;
        for (int i = 0; i < gap; i++)
            rat_mul(result, result, x);
        rat_add(result, result, &t->coeff);
        prev_exp = t->exponent;
        t = t->next;
    }
    for (int i = 0; i < prev_exp; i++)
        rat_mul(result, result, x);
}

// collect all positive integer divisors of |n|
static long* divisors_long(long n, int* count) {
    if (n < 0) n = -n;
    if (n == 0) { *count = 0; return NULL; }
    int cap = 16;
    long* divs = malloc(cap * sizeof(long));
    *count = 0;
    for (long i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            if (*count + 2 > cap) { cap *= 2; divs = realloc(divs, cap * sizeof(long)); }
            divs[(*count)++] = i;
            if (i != n / i)
                divs[(*count)++] = n / i;
        }
    }
    return divs;
}

// clear denominators to get integer coefficients; returns lcm of denominators
static Polynomial* clear_denominators(const Polynomial* poly, mpz_t lcm_den) {
    mpz_set_si(lcm_den, 1);
    Term* t = poly->head;
    mpz_t g;
    mpz_init(g);
    while (t) {
        mpz_gcd(g, lcm_den, mpq_denref(t->coeff.val));
        mpz_divexact(lcm_den, lcm_den, g);
        mpz_mul(lcm_den, lcm_den, mpq_denref(t->coeff.val));
        t = t->next;
    }
    Polynomial* result = create_polynomial();
    t = poly->head;
    mpz_t new_num, d;
    mpz_init(new_num);
    mpz_init(d);
    while (t) {
        mpz_divexact(d, lcm_den, mpq_denref(t->coeff.val));
        mpz_mul(new_num, mpq_numref(t->coeff.val), d);
        Rational r;
        rat_init(&r);
        mpq_set_z(r.val, new_num);
        add_term(result, &r, t->exponent);
        rat_clear(&r);
        t = t->next;
    }
    mpz_clear(new_num);
    mpz_clear(d);
    mpz_clear(g);
    return result;
}

static Rational* find_rational_roots(const Polynomial* poly, int* num_roots) {
    *num_roots = 0;
    if (poly->head == NULL) return NULL;

    mpz_t lcm_den;
    mpz_init(lcm_den);
    Polynomial* int_poly = clear_denominators(poly, lcm_den);
    mpz_clear(lcm_den);

    long leading = rat_num_si(&int_poly->head->coeff);
    long constant = 0;
    Term* t = int_poly->head;
    while (t) {
        if (t->exponent == 0) { constant = rat_num_si(&t->coeff); break; }
        t = t->next;
    }

    if (constant == 0) {
        int cap = 8;
        Rational* roots = malloc(cap * sizeof(Rational));
        roots[(*num_roots)++] = rat_zero_val();
        Rational zero_r = rat_zero_val();
        Polynomial* reduced = synthetic_divide(int_poly, &zero_r);
        rat_clear(&zero_r);
        int more_count;
        Rational* more_roots = find_rational_roots(reduced, &more_count);
        for (int i = 0; i < more_count; i++) {
            if (*num_roots >= cap) { cap *= 2; roots = realloc(roots, cap * sizeof(Rational)); }
            roots[*num_roots] = more_roots[i]; // transfer ownership
            (*num_roots)++;
        }
        free(more_roots);
        free_polynomial(reduced);
        free_polynomial(int_poly);
        return roots;
    }

    int num_p, num_q;
    long* p_divs = divisors_long(constant, &num_p);
    long* q_divs = divisors_long(leading, &num_q);

    int cap = 16;
    Rational* roots = malloc(cap * sizeof(Rational));
    Polynomial* current = poly_copy(int_poly);

    for (int i = 0; i < num_p && !poly_is_zero(current); i++) {
        for (int j = 0; j < num_q && !poly_is_zero(current); j++) {
            Rational candidates[2];
            candidates[0] = rat_from_si(p_divs[i], q_divs[j]);
            candidates[1] = rat_from_si(-p_divs[i], q_divs[j]);
            for (int c = 0; c < 2 && !poly_is_zero(current); c++) {
                Rational val;
                rat_init(&val);
                poly_eval_rational(&val, current, &candidates[c]);
                if (rat_is_zero(&val)) {
                    if (*num_roots >= cap) { cap *= 2; roots = realloc(roots, cap * sizeof(Rational)); }
                    rat_init_set(&roots[*num_roots], &candidates[c]);
                    (*num_roots)++;
                    Polynomial* reduced = synthetic_divide(current, &candidates[c]);
                    free_polynomial(current);
                    current = reduced;
                    c--;
                }
                rat_clear(&val);
            }
            rat_clear(&candidates[0]);
            rat_clear(&candidates[1]);
        }
    }

    free(p_divs);
    free(q_divs);
    free_polynomial(current);
    free_polynomial(int_poly);
    return roots;
}

Factorization* factorize(const Polynomial* poly) {
    Factorization* f = malloc(sizeof(Factorization));
    f->factors = NULL;
    f->multiplicities = NULL;
    f->count = 0;
    f->content = rat_one_val();

    if (poly->head == NULL) {
        rat_clear(&f->content);
        f->content = rat_zero_val();
        return f;
    }

    rat_clear(&f->content);
    f->content = poly_content(poly);
    if (rat_is_negative(&poly->head->coeff)) {
        rat_neg(&f->content, &f->content);
    }

    Polynomial* primitive = create_polynomial();
    Term* t = poly->head;
    while (t) {
        Rational div_r;
        rat_init(&div_r);
        rat_div(&div_r, &t->coeff, &f->content);
        add_term(primitive, &div_r, t->exponent);
        rat_clear(&div_r);
        t = t->next;
    }

    int num_roots;
    Rational* roots = find_rational_roots(primitive, &num_roots);

    int cap = 8;
    f->factors = malloc(cap * sizeof(Polynomial*));
    f->multiplicities = malloc(cap * sizeof(int));

    Polynomial* remainder = poly_copy(primitive);

    for (int i = 0; i < num_roots; i++) {
        bool found = false;
        for (int j = 0; j < f->count; j++) {
            Term* ft = f->factors[j]->head;
            if (ft && ft->exponent == 1 && ft->next) {
                Rational neg_root;
                rat_init(&neg_root);
                rat_neg(&neg_root, &roots[i]);
                if (rat_eq(&ft->next->coeff, &neg_root)) {
                    f->multiplicities[j]++;
                    found = true;
                }
                rat_clear(&neg_root);
                if (found) break;
            }
        }
        if (!found) {
            if (f->count >= cap) {
                cap *= 2;
                f->factors = realloc(f->factors, cap * sizeof(Polynomial*));
                f->multiplicities = realloc(f->multiplicities, cap * sizeof(int));
            }
            Polynomial* factor = create_polynomial();
            Rational one_r = rat_one_val();
            add_term(factor, &one_r, 1);
            rat_clear(&one_r);
            Rational neg_root;
            rat_init(&neg_root);
            rat_neg(&neg_root, &roots[i]);
            add_term(factor, &neg_root, 0);
            rat_clear(&neg_root);
            f->factors[f->count] = factor;
            f->multiplicities[f->count] = 1;
            f->count++;
        }

        Polynomial* reduced = synthetic_divide(remainder, &roots[i]);
        free_polynomial(remainder);
        remainder = reduced;
    }

    if (!poly_is_zero(remainder) && polynomial_degree(remainder) > 0) {
        if (f->count >= cap) {
            cap *= 2;
            f->factors = realloc(f->factors, cap * sizeof(Polynomial*));
            f->multiplicities = realloc(f->multiplicities, cap * sizeof(int));
        }
        f->factors[f->count] = poly_copy(remainder);
        f->multiplicities[f->count] = 1;
        f->count++;
    }

    for (int i = 0; i < num_roots; i++) rat_clear(&roots[i]);
    free(roots);
    free_polynomial(remainder);
    free_polynomial(primitive);
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
    size_t cap = 128, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';

    Rational neg_one = rat_from_int(-1);
    if (!rat_is_one(&f->content) && !rat_eq(&f->content, &neg_one)) {
        char* cs = rat_to_string(&f->content);
        while (len + strlen(cs) + 2 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, cs);
        len += strlen(cs);
        free(cs);
    } else if (rat_eq(&f->content, &neg_one)) {
        strcat(buf, "-");
        len += 1;
    }
    rat_clear(&neg_one);

    if (f->count == 0) {
        if (len == 0) {
            char* cs = rat_to_string(&f->content);
            while (len + strlen(cs) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
            strcat(buf, cs);
            free(cs);
        }
        return buf;
    }

    for (int i = 0; i < f->count; i++) {
        char* ps = poly_to_string(f->factors[i]);
        char factor_buf[256];
        if (f->multiplicities[i] > 1)
            snprintf(factor_buf, sizeof(factor_buf), "(%s)^%d", ps, f->multiplicities[i]);
        else
            snprintf(factor_buf, sizeof(factor_buf), "(%s)", ps);
        free(ps);
        size_t flen = strlen(factor_buf);
        while (len + flen + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, factor_buf);
        len += flen;
    }

    return buf;
}

Expr* factorization_to_expr(const Factorization* f) {
    (void)f;
    return NULL;
}
