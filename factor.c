#include "factor.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// GCD of all coefficients in the polynomial
// gcd(a/b, c/d) = gcd(a,c) / lcm(b,d)
Rational poly_content(const Polynomial* poly) {
    if (poly->head == NULL) return rat_zero();
    Rational content = rat_abs(poly->head->coeff);
    Term* t = poly->head->next;
    while (t) {
        Rational abs_c = rat_abs(t->coeff);
        rat_int_t num_gcd = gcd_int(content.num, abs_c.num);
        rat_int_t den_lcm = content.den / gcd_int(content.den, abs_c.den) * abs_c.den;
        content = rat_create(num_gcd, den_lcm);
        t = t->next;
    }
    return content;
}

// divide all coefficients by content, yielding a primitive polynomial
Polynomial* poly_primitive_part(const Polynomial* poly) {
    Rational c = poly_content(poly);
    if (rat_is_zero(c)) return poly_copy(poly);
    Polynomial* result = create_polynomial();
    Term* t = poly->head;
    while (t) {
        add_term(result, rat_div(t->coeff, c), t->exponent);
        t = t->next;
    }
    return result;
}

// divide poly by (x - root) using synthetic division
Polynomial* synthetic_divide(const Polynomial* poly, Rational root) {
    int deg = polynomial_degree(poly);
    if (deg <= 0) return create_polynomial();

    // convert sparse linked list to dense coefficient array (descending)
    Rational* coeffs = calloc(deg + 1, sizeof(Rational));
    for (int i = 0; i <= deg; i++)
        coeffs[i] = rat_zero();
    Term* t = poly->head;
    while (t) {
        coeffs[deg - t->exponent] = t->coeff;
        t = t->next;
    }

    // synthetic division
    Rational* result = calloc(deg, sizeof(Rational));
    result[0] = coeffs[0];
    for (int i = 1; i < deg; i++)
        result[i] = rat_add(coeffs[i], rat_mul(result[i - 1], root));

    Polynomial* quotient = create_polynomial();
    for (int i = 0; i < deg; i++) {
        if (!rat_is_zero(result[i]))
            add_term(quotient, result[i], deg - 1 - i);
    }

    free(coeffs);
    free(result);
    return quotient;
}

// evaluate polynomial at a rational point (exact)
static Rational poly_eval_rational(const Polynomial* poly, Rational x) {
    if (poly->head == NULL) return rat_zero();
    Term* t = poly->head;
    Rational result = t->coeff;
    int prev_exp = t->exponent;
    t = t->next;
    while (t) {
        int gap = prev_exp - t->exponent;
        for (int i = 0; i < gap; i++)
            result = rat_mul(result, x);
        result = rat_add(result, t->coeff);
        prev_exp = t->exponent;
        t = t->next;
    }
    for (int i = 0; i < prev_exp; i++)
        result = rat_mul(result, x);
    return result;
}

// collect all integer divisors of |n| (including 1 and |n|)
static rat_int_t* divisors(rat_int_t n, int* count) {
    if (n < 0) n = -n;
    if (n == 0) { *count = 0; return NULL; }
    int cap = 16;
    rat_int_t* divs = malloc(cap * sizeof(rat_int_t));
    *count = 0;
    for (rat_int_t i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            if (*count + 2 > cap) { cap *= 2; divs = realloc(divs, cap * sizeof(rat_int_t)); }
            divs[(*count)++] = i;
            if (i != n / i)
                divs[(*count)++] = n / i;
        }
    }
    return divs;
}

// clear denominators: multiply poly so all coefficients are integers
// returns the cleared polynomial (caller frees) and sets *lcm_den
static Polynomial* clear_denominators(const Polynomial* poly, rat_int_t* lcm_den) {
    *lcm_den = 1;
    Term* t = poly->head;
    while (t) {
        rat_int_t g = gcd_int(*lcm_den, t->coeff.den);
        *lcm_den = (*lcm_den / g) * t->coeff.den;
        t = t->next;
    }
    Polynomial* result = create_polynomial();
    t = poly->head;
    while (t) {
        rat_int_t new_num = t->coeff.num * (*lcm_den / t->coeff.den);
        add_term(result, rat_from_int(new_num), t->exponent);
        t = t->next;
    }
    return result;
}

// find all rational roots using the rational root theorem
static Rational* find_rational_roots(const Polynomial* poly, int* num_roots) {
    *num_roots = 0;
    if (poly->head == NULL) return NULL;

    // clear denominators to get integer coefficients
    rat_int_t lcm_den;
    Polynomial* int_poly = clear_denominators(poly, &lcm_den);

    // find leading and trailing coefficients
    rat_int_t leading = int_poly->head->coeff.num;
    // find constant term
    rat_int_t constant = 0;
    Term* t = int_poly->head;
    while (t) {
        if (t->exponent == 0) { constant = t->coeff.num; break; }
        t = t->next;
    }

    if (constant == 0) {
        // x = 0 is a root
        int cap = 8;
        Rational* roots = malloc(cap * sizeof(Rational));
        roots[(*num_roots)++] = rat_zero();
        // divide out x and find more roots
        Polynomial* reduced = synthetic_divide(int_poly, rat_zero());
        int more_count;
        Rational* more_roots = find_rational_roots(reduced, &more_count);
        for (int i = 0; i < more_count; i++) {
            if (*num_roots >= cap) { cap *= 2; roots = realloc(roots, cap * sizeof(Rational)); }
            roots[(*num_roots)++] = more_roots[i];
        }
        free(more_roots);
        free_polynomial(reduced);
        free_polynomial(int_poly);
        return roots;
    }

    int num_p, num_q;
    rat_int_t* p_divs = divisors(constant, &num_p);
    rat_int_t* q_divs = divisors(leading, &num_q);

    int cap = 16;
    Rational* roots = malloc(cap * sizeof(Rational));
    Polynomial* current = poly_copy(int_poly);

    // test all p/q and -p/q candidates
    for (int i = 0; i < num_p && !poly_is_zero(current); i++) {
        for (int j = 0; j < num_q && !poly_is_zero(current); j++) {
            Rational candidates[2] = {
                rat_create(p_divs[i], q_divs[j]),
                rat_create(-p_divs[i], q_divs[j])
            };
            for (int c = 0; c < 2 && !poly_is_zero(current); c++) {
                Rational val = poly_eval_rational(current, candidates[c]);
                if (rat_is_zero(val)) {
                    if (*num_roots >= cap) { cap *= 2; roots = realloc(roots, cap * sizeof(Rational)); }
                    roots[(*num_roots)++] = candidates[c];
                    // divide out the root and check for repeated roots
                    Polynomial* reduced = synthetic_divide(current, candidates[c]);
                    free_polynomial(current);
                    current = reduced;
                    // check if same root appears again
                    c--;
                }
            }
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
    f->content = rat_one();

    if (poly->head == NULL) {
        f->content = rat_zero();
        return f;
    }

    // extract content
    f->content = poly_content(poly);
    if (rat_is_negative(poly->head->coeff))
        f->content = rat_neg(f->content);

    Polynomial* primitive = create_polynomial();
    Term* t = poly->head;
    while (t) {
        add_term(primitive, rat_div(t->coeff, f->content), t->exponent);
        t = t->next;
    }

    // find rational roots
    int num_roots;
    Rational* roots = find_rational_roots(primitive, &num_roots);

    // build factor list
    int cap = 8;
    f->factors = malloc(cap * sizeof(Polynomial*));
    f->multiplicities = malloc(cap * sizeof(int));

    Polynomial* remainder = poly_copy(primitive);

    for (int i = 0; i < num_roots; i++) {
        // check if this root is already in the factor list
        bool found = false;
        for (int j = 0; j < f->count; j++) {
            // compare: factor is (x - root), which is degree 1 with coeff 1 and constant -root
            Term* ft = f->factors[j]->head;
            if (ft && ft->exponent == 1 && ft->next &&
                rat_eq(ft->next->coeff, rat_neg(roots[i]))) {
                f->multiplicities[j]++;
                found = true;
                break;
            }
        }
        if (!found) {
            if (f->count >= cap) {
                cap *= 2;
                f->factors = realloc(f->factors, cap * sizeof(Polynomial*));
                f->multiplicities = realloc(f->multiplicities, cap * sizeof(int));
            }
            // create factor (x - root)
            Polynomial* factor = create_polynomial();
            add_term(factor, rat_one(), 1);
            add_term(factor, rat_neg(roots[i]), 0);
            f->factors[f->count] = factor;
            f->multiplicities[f->count] = 1;
            f->count++;
        }

        Polynomial* reduced = synthetic_divide(remainder, roots[i]);
        free_polynomial(remainder);
        remainder = reduced;
    }

    // if remainder has degree > 0, add it as an irreducible factor
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
    free(f);
}

char* factorization_to_string(const Factorization* f) {
    size_t cap = 128, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';

    // content prefix
    if (!rat_is_one(f->content) && !rat_eq(f->content, rat_from_int(-1))) {
        char* cs = rat_to_string(f->content);
        while (len + strlen(cs) + 2 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, cs);
        len += strlen(cs);
        free(cs);
    } else if (rat_eq(f->content, rat_from_int(-1))) {
        strcat(buf, "-");
        len += 1;
    }

    if (f->count == 0) {
        if (len == 0) {
            char* cs = rat_to_string(f->content);
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
    return NULL; // will be used by CLI
}
