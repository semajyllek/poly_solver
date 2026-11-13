#include "polynomial.h"
#include "strbuf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Term* create_term(const Rational* coeff, int exponent) {
    Term* term = malloc(sizeof(Term));
    rat_init_set(&term->coeff, coeff);
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

void add_term(Polynomial* poly, const Rational* coeff, int exponent) {
    if (rat_is_zero(coeff)) return;

    Term* new_term = create_term(coeff, exponent);
    if (poly->head == NULL) {
        poly->head = new_term;
        return;
    }

    Term* current = poly->head;
    Term* previous = NULL;

    while (current != NULL && current->exponent > exponent) {
        previous = current;
        current = current->next;
    }

    if (current != NULL && current->exponent == exponent) {
        rat_add(&current->coeff, &current->coeff, coeff);
        rat_clear(&new_term->coeff);
        free(new_term);
        if (rat_is_zero(&current->coeff)) {
            if (previous == NULL)
                poly->head = current->next;
            else
                previous->next = current->next;
            rat_clear(&current->coeff);
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
        add_term(copy, &current->coeff, current->exponent);
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
        add_term(result, &current->coeff, current->exponent);
        current = current->next;
    }
    current = poly2->head;
    while (current != NULL) {
        add_term(result, &current->coeff, current->exponent);
        current = current->next;
    }
    return result;
}

Polynomial* poly_subtract(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* result = create_polynomial();
    Term* current = poly1->head;
    while (current != NULL) {
        add_term(result, &current->coeff, current->exponent);
        current = current->next;
    }
    current = poly2->head;
    while (current != NULL) {
        Rational neg;
        rat_init(&neg);
        rat_neg(&neg, &current->coeff);
        add_term(result, &neg, current->exponent);
        rat_clear(&neg);
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
            Rational prod;
            rat_init(&prod);
            rat_mul(&prod, &t1->coeff, &t2->coeff);
            add_term(result, &prod, t1->exponent + t2->exponent);
            rat_clear(&prod);
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
        Rational coeff;
        rat_init(&coeff);
        rat_div(&coeff, &remainder->head->coeff, &divisor->head->coeff);
        int exp = remainder->head->exponent - divisor->head->exponent;

        add_term(result.quotient, &coeff, exp);

        Polynomial* term_poly = create_polynomial();
        add_term(term_poly, &coeff, exp);
        Polynomial* product = poly_multiply(term_poly, divisor);
        Polynomial* new_remainder = poly_subtract(remainder, product);

        free_polynomial(remainder);
        free_polynomial(product);
        free_polynomial(term_poly);
        rat_clear(&coeff);
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
    if (rat_is_one(&poly->head->coeff)) return;
    Rational lc;
    rat_init_set(&lc, &poly->head->coeff);
    Term* current = poly->head;
    while (current != NULL) {
        Rational tmp;
        rat_init(&tmp);
        rat_div(&tmp, &current->coeff, &lc);
        rat_set(&current->coeff, &tmp);
        rat_clear(&tmp);
        current = current->next;
    }
    rat_clear(&lc);
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

Polynomial* poly_constant(const Rational* c) {
    Polynomial* p = create_polynomial();
    add_term(p, c, 0);
    return p;
}

Polynomial* poly_pow(const Polynomial* base, unsigned exp) {
    Rational one = rat_one_val();
    Polynomial* result = poly_constant(&one);
    rat_clear(&one);
    Polynomial* sq = poly_copy(base);
    while (exp > 0) {
        if (exp & 1u) { // fold the current square into the result
            Polynomial* t = poly_multiply(result, sq);
            free_polynomial(result);
            result = t;
        }
        exp >>= 1;
        if (exp) {
            Polynomial* t = poly_multiply(sq, sq);
            free_polynomial(sq);
            sq = t;
        }
    }
    free_polynomial(sq);
    return result;
}

Polynomial* poly_scale(const Polynomial* p, const Rational* c) {
    Polynomial* r = create_polynomial();
    Rational t;
    rat_init(&t);
    for (Term* term = p->head; term; term = term->next) {
        rat_mul(&t, &term->coeff, c);
        add_term(r, &t, term->exponent);
    }
    rat_clear(&t);
    return r;
}

// next = old - q * cur (one row of the extended Euclidean recurrence)
static Polynomial* euclid_row(const Polynomial* old, const Polynomial* q, const Polynomial* cur) {
    Polynomial* qc = poly_multiply(q, cur);
    Polynomial* r = poly_subtract(old, qc);
    free_polynomial(qc);
    return r;
}

PolyBezout poly_xgcd(const Polynomial* a, const Polynomial* b) {
    Rational one = rat_one_val();
    Polynomial *r0 = poly_copy(a), *r1 = poly_copy(b);
    Polynomial *s0 = poly_constant(&one), *s1 = create_polynomial();
    Polynomial *t0 = create_polynomial(), *t1 = poly_constant(&one);
    while (!poly_is_zero(r1)) {
        PolyDivResult dr = poly_divmod(r0, r1);
        Polynomial* s2 = euclid_row(s0, dr.quotient, s1);
        Polynomial* t2 = euclid_row(t0, dr.quotient, t1);
        free_polynomial(dr.quotient); free_polynomial(r0); free_polynomial(s0); free_polynomial(t0);
        r0 = r1; r1 = dr.remainder;
        s0 = s1; s1 = s2;
        t0 = t1; t1 = t2;
    }
    // normalize so the gcd is monic
    Rational inv;
    rat_init(&inv);
    if (r0->head) mpq_inv(inv.val, r0->head->coeff.val);
    else mpq_set_si(inv.val, 1, 1);
    PolyBezout bz = { poly_scale(s0, &inv), poly_scale(t0, &inv), poly_scale(r0, &inv) };
    free_polynomial(r0); free_polynomial(r1); free_polynomial(s0);
    free_polynomial(s1); free_polynomial(t0); free_polynomial(t1);
    rat_clear(&inv); rat_clear(&one);
    return bz;
}

void free_bezout(PolyBezout* bz) {
    free_polynomial(bz->s);
    free_polynomial(bz->t);
    free_polynomial(bz->g);
}

void poly_coeff(Rational* out, const Polynomial* poly, int degree) {
    for (Term* t = poly->head; t; t = t->next) {
        if (t->exponent == degree) { rat_set(out, &t->coeff); return; }
        if (t->exponent < degree) break; // terms are sorted descending
    }
    mpq_set_si(out->val, 0, 1);
}

void poly_eval_rat(Rational* out, const Polynomial* poly, const Rational* x) {
    mpq_set_si(out->val, 0, 1);
    int prev = poly->head ? poly->head->exponent : 0;
    for (Term* t = poly->head; t; t = t->next) {
        for (int i = t->exponent; i < prev; i++) rat_mul(out, out, x); // skip the gap
        rat_add(out, out, &t->coeff);
        prev = t->exponent;
    }
    for (int i = 0; i < prev; i++) rat_mul(out, out, x);
}

double poly_evaluate(const Polynomial* poly, double x) {
    if (poly->head == NULL) return 0.0;

    Term* current = poly->head;
    double result = rat_to_double(&current->coeff);
    int prev_exp = current->exponent;
    current = current->next;

    while (current != NULL) {
        int gap = prev_exp - current->exponent;
        for (int i = 0; i < gap; i++)
            result *= x;
        result += rat_to_double(&current->coeff);
        prev_exp = current->exponent;
        current = current->next;
    }

    for (int i = 0; i < prev_exp; i++)
        result *= x;

    return result;
}

Polynomial* poly_derivative(const Polynomial* poly) {
    Polynomial* derivative = create_polynomial();
    Term* current = poly->head;
    while (current != NULL) {
        if (current->exponent > 0) {
            Rational exp_rat = rat_from_int(current->exponent);
            Rational new_coeff;
            rat_init(&new_coeff);
            rat_mul(&new_coeff, &current->coeff, &exp_rat);
            add_term(derivative, &new_coeff, current->exponent - 1);
            rat_clear(&new_coeff);
            rat_clear(&exp_rat);
        }
        current = current->next;
    }
    return derivative;
}

Polynomial* poly_integral(const Polynomial* poly) {
    Polynomial* integral = create_polynomial();
    Term* current = poly->head;
    while (current != NULL) {
        Rational divisor = rat_from_int(current->exponent + 1);
        Rational new_coeff;
        rat_init(&new_coeff);
        rat_div(&new_coeff, &current->coeff, &divisor);
        add_term(integral, &new_coeff, current->exponent + 1);
        rat_clear(&new_coeff);
        rat_clear(&divisor);
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
        rat_clear(&temp->coeff);
        free(temp);
    }
    free(poly);
}

void print_polynomial(const Polynomial* poly) {
    char* s = poly_to_string(poly);
    printf("%s", s);
    free(s);
}

typedef struct {
    const Rational* abs_c; // |coefficient|
    int e;                 // exponent
    const char* var;       // variable name
} TermFmt;

// append one term (sign handled by caller) as "c", "x^e", or "cx^e"
static void append_term_body(StrBuf* sb, const TermFmt* t) {
    if (t->e == 0) { sb_append_owned(sb, rat_to_string(t->abs_c)); return; }
    if (!rat_is_one(t->abs_c)) sb_append_owned(sb, rat_to_string(t->abs_c));
    sb_append(sb, t->var);
    if (t->e != 1) sb_appendf(sb, "^%d", t->e);
}

char* poly_to_string_var(const Polynomial* poly, const char* var) {
    StrBuf sb = sb_new();
    if (poly->head == NULL) { sb_append(&sb, "0"); return sb_finish(&sb); }
    Rational c;
    rat_init(&c);
    for (Term* t = poly->head; t; t = t->next) {
        bool neg = rat_is_negative(&t->coeff);
        if (t != poly->head) sb_append(&sb, neg ? " - " : " + ");
        else if (neg) sb_append(&sb, "-");
        rat_abs(&c, &t->coeff);
        TermFmt fmt = { &c, t->exponent, var };
        append_term_body(&sb, &fmt);
    }
    rat_clear(&c);
    return sb_finish(&sb);
}

char* poly_to_string(const Polynomial* poly) {
    return poly_to_string_var(poly, "x");
}

QuadRoots poly_quad_roots(const Polynomial* f) {
    QuadRoots q;
    Rational a, b, c, d, t;
    rat_init(&q.p); rat_init(&q.r); mpz_init(q.m);
    rat_init(&a); rat_init(&b); rat_init(&c); rat_init(&d); rat_init(&t);
    poly_coeff(&a, f, 2); poly_coeff(&b, f, 1); poly_coeff(&c, f, 0);
    rat_mul(&d, &b, &b);                       // D = b^2 - 4ac
    rat_mul(&t, &a, &c);
    mpq_mul_2exp(t.val, t.val, 2);
    rat_sub(&d, &d, &t);
    q.is_complex = rat_is_negative(&d);
    rat_abs(&d, &d);
    rat_mul(&t, &a, &a);                       // |D| / (4a^2)
    mpq_mul_2exp(t.val, t.val, 2);
    rat_div(&d, &d, &t);
    rat_sqrt_split(&d, &q.r, q.m);
    rat_div(&q.p, &b, &a);                     // p = -b / 2a
    mpq_div_2exp(q.p.val, q.p.val, 1);
    rat_neg(&q.p, &q.p);
    rat_clear(&a); rat_clear(&b); rat_clear(&c); rat_clear(&d); rat_clear(&t);
    return q;
}

void quad_roots_clear(QuadRoots* q) {
    rat_clear(&q->p); rat_clear(&q->r); mpz_clear(q->m);
}
