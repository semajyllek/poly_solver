#include "polynomial.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Term* create_term(Rational coeff, int exponent) {
    Term* term = malloc(sizeof(Term));
    term->coeff = coeff;
    term->exponent = exponent;
    term->var_id = 0;
    term->next = NULL;
    return term;
}

Polynomial* create_polynomial(void) {
    Polynomial* poly = malloc(sizeof(Polynomial));
    poly->head = NULL;
    return poly;
}

void add_term(Polynomial* poly, Rational coeff, int exponent) {
    if (rat_is_zero(coeff)) return;

    Term* new_term = create_term(coeff, exponent);
    if (poly->head == NULL) {
        poly->head = new_term;
        return;
    }

    Term* current = poly->head;
    Term* previous = NULL;

    // find correct position (descending exponent order)
    while (current != NULL && current->exponent > exponent) {
        previous = current;
        current = current->next;
    }

    if (current != NULL && current->exponent == exponent) {
        // combine like terms
        current->coeff = rat_add(current->coeff, coeff);
        free(new_term);
        if (rat_is_zero(current->coeff)) {
            if (previous == NULL)
                poly->head = current->next;
            else
                previous->next = current->next;
            free(current);
        }
    } else {
        new_term->next = current;
        if (previous == NULL)
            poly->head = new_term;
        else
            previous->next = new_term;
    }
}

Polynomial* poly_copy(const Polynomial* poly) {
    Polynomial* copy = create_polynomial();
    Term* current = poly->head;
    while (current != NULL) {
        add_term(copy, current->coeff, current->exponent);
        current = current->next;
    }
    return copy;
}

bool poly_is_zero(const Polynomial* poly) {
    return poly->head == NULL;
}

int polynomial_degree(const Polynomial* poly) {
    if (poly->head == NULL) return -1;
    return poly->head->exponent;
}

Polynomial* poly_add(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* result = create_polynomial();
    Term* current = poly1->head;
    while (current != NULL) {
        add_term(result, current->coeff, current->exponent);
        current = current->next;
    }
    current = poly2->head;
    while (current != NULL) {
        add_term(result, current->coeff, current->exponent);
        current = current->next;
    }
    return result;
}

Polynomial* poly_subtract(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* result = create_polynomial();
    Term* current = poly1->head;
    while (current != NULL) {
        add_term(result, current->coeff, current->exponent);
        current = current->next;
    }
    current = poly2->head;
    while (current != NULL) {
        add_term(result, rat_neg(current->coeff), current->exponent);
        current = current->next;
    }
    return result;
}

Polynomial* poly_multiply(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* result = create_polynomial();
    Term* t1 = poly1->head;
    while (t1 != NULL) {
        Term* t2 = poly2->head;
        while (t2 != NULL) {
            add_term(result, rat_mul(t1->coeff, t2->coeff), t1->exponent + t2->exponent);
            t2 = t2->next;
        }
        t1 = t1->next;
    }
    return result;
}

PolyDivResult poly_divmod(const Polynomial* dividend, const Polynomial* divisor) {
    PolyDivResult result = {create_polynomial(), NULL};

    if (divisor == NULL || divisor->head == NULL) {
        fprintf(stderr, "Error: division by zero polynomial.\n");
        result.remainder = create_polynomial();
        return result;
    }

    Polynomial* remainder = poly_copy(dividend);

    while (!poly_is_zero(remainder) &&
           remainder->head->exponent >= divisor->head->exponent) {
        Rational coeff = rat_div(remainder->head->coeff, divisor->head->coeff);
        int exp = remainder->head->exponent - divisor->head->exponent;

        add_term(result.quotient, coeff, exp);

        // subtract coeff * x^exp * divisor from remainder
        Polynomial* term_poly = create_polynomial();
        add_term(term_poly, coeff, exp);
        Polynomial* product = poly_multiply(term_poly, divisor);
        Polynomial* new_remainder = poly_subtract(remainder, product);

        free_polynomial(remainder);
        free_polynomial(product);
        free_polynomial(term_poly);
        remainder = new_remainder;
    }

    result.remainder = remainder;
    return result;
}

Polynomial* poly_divide(const Polynomial* dividend, const Polynomial* divisor) {
    PolyDivResult dr = poly_divmod(dividend, divisor);
    free_polynomial(dr.remainder);
    return dr.quotient;
}

Polynomial* poly_mod(const Polynomial* dividend, const Polynomial* divisor) {
    PolyDivResult dr = poly_divmod(dividend, divisor);
    free_polynomial(dr.quotient);
    return dr.remainder;
}

void poly_make_monic(Polynomial* poly) {
    if (poly->head == NULL) return;
    Rational lc = poly->head->coeff;
    if (rat_is_one(lc)) return;
    Term* current = poly->head;
    while (current != NULL) {
        current->coeff = rat_div(current->coeff, lc);
        current = current->next;
    }
}

Polynomial* poly_gcd(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* r0 = poly_copy(poly1);
    Polynomial* r1 = poly_copy(poly2);

    while (!poly_is_zero(r1)) {
        Polynomial* r2 = poly_mod(r0, r1);
        free_polynomial(r0);
        r0 = r1;
        r1 = r2;
    }

    free_polynomial(r1);
    poly_make_monic(r0);
    return r0;
}

Polynomial* poly_lcm(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* product = poly_multiply(poly1, poly2);
    Polynomial* gcd = poly_gcd(poly1, poly2);
    Polynomial* lcm = poly_divide(product, gcd);
    free_polynomial(product);
    free_polynomial(gcd);
    return lcm;
}

// Horner's method for evaluation
double poly_evaluate(const Polynomial* poly, double x) {
    if (poly->head == NULL) return 0.0;

    Term* current = poly->head;
    double result = rat_to_double(current->coeff);
    int prev_exp = current->exponent;
    current = current->next;

    while (current != NULL) {
        // multiply by x for each degree gap
        int gap = prev_exp - current->exponent;
        for (int i = 0; i < gap; i++)
            result *= x;
        result += rat_to_double(current->coeff);
        prev_exp = current->exponent;
        current = current->next;
    }

    // multiply by x for remaining degrees down to 0
    for (int i = 0; i < prev_exp; i++)
        result *= x;

    return result;
}

Polynomial* poly_derivative(const Polynomial* poly) {
    Polynomial* derivative = create_polynomial();
    Term* current = poly->head;
    while (current != NULL) {
        if (current->exponent > 0) {
            Rational new_coeff = rat_mul(current->coeff, rat_from_int(current->exponent));
            add_term(derivative, new_coeff, current->exponent - 1);
        }
        current = current->next;
    }
    return derivative;
}

Polynomial* poly_integral(const Polynomial* poly) {
    Polynomial* integral = create_polynomial();
    Term* current = poly->head;
    while (current != NULL) {
        Rational new_coeff = rat_div(current->coeff, rat_from_int(current->exponent + 1));
        add_term(integral, new_coeff, current->exponent + 1);
        current = current->next;
    }
    return integral;
}

double poly_definite_integral(const Polynomial* poly, double lower, double upper) {
    Polynomial* integral = poly_integral(poly);
    double result = poly_evaluate(integral, upper) - poly_evaluate(integral, lower);
    free_polynomial(integral);
    return result;
}

void free_polynomial(Polynomial* poly) {
    if (!poly) return;
    Term* current = poly->head;
    while (current != NULL) {
        Term* temp = current;
        current = current->next;
        free(temp);
    }
    free(poly);
}

void print_polynomial(const Polynomial* poly) {
    char* s = poly_to_string(poly);
    printf("%s", s);
    free(s);
}

char* poly_to_string(const Polynomial* poly) {
    if (poly->head == NULL) {
        char* s = malloc(2);
        strcpy(s, "0");
        return s;
    }

    // build string incrementally
    size_t cap = 128, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';

    Term* current = poly->head;
    bool first = true;

    while (current != NULL) {
        char term_buf[128];
        Rational c = current->coeff;
        int e = current->exponent;

        if (!first) {
            if (rat_is_negative(c)) {
                strcat(buf, " - ");
                len += 3;
                c = rat_neg(c);
            } else {
                strcat(buf, " + ");
                len += 3;
            }
        } else if (rat_is_negative(c)) {
            strcat(buf, "-");
            len += 1;
            c = rat_neg(c);
        }

        bool coeff_is_one = rat_is_one(c);

        if (e == 0) {
            // constant term — always show coefficient
            char* cs = rat_to_string(c);
            snprintf(term_buf, sizeof(term_buf), "%s", cs);
            free(cs);
        } else if (coeff_is_one) {
            if (e == 1)
                snprintf(term_buf, sizeof(term_buf), "x");
            else
                snprintf(term_buf, sizeof(term_buf), "x^%d", e);
        } else {
            char* cs = rat_to_string(c);
            if (e == 1)
                snprintf(term_buf, sizeof(term_buf), "%sx", cs);
            else
                snprintf(term_buf, sizeof(term_buf), "%sx^%d", cs, e);
            free(cs);
        }

        size_t tlen = strlen(term_buf);
        while (len + tlen + 1 > cap) {
            cap *= 2;
            buf = realloc(buf, cap);
        }
        strcat(buf, term_buf);
        len += tlen;

        first = false;
        current = current->next;
    }

    return buf;
}
