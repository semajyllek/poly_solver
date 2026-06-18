#include "solver.h"
#include "factor.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static bool is_perfect_square_rat(const Rational* r, Rational* root) {
    if (rat_is_negative(r)) return false;
    if (rat_is_zero(r)) { mpq_set_si(root->val, 0, 1); return true; }
    // check if both num and den are perfect squares
    if (!mpz_perfect_square_p(mpq_numref(r->val))) return false;
    if (!mpz_perfect_square_p(mpq_denref(r->val))) return false;
    mpz_t sn, sd;
    mpz_init(sn); mpz_init(sd);
    mpz_sqrt(sn, mpq_numref(r->val));
    mpz_sqrt(sd, mpq_denref(r->val));
    mpq_set_num(root->val, sn);
    mpq_set_den(root->val, sd);
    mpq_canonicalize(root->val);
    mpz_clear(sn); mpz_clear(sd);
    return true;
}

static void solve_quadratic(const Rational* a, const Rational* b, const Rational* c, Solutions* sol) {
    // D = b^2 - 4ac
    Rational b2, four, ac, four_ac, d;
    rat_init(&b2); rat_init(&four); rat_init(&ac); rat_init(&four_ac); rat_init(&d);
    rat_mul(&b2, b, b);
    mpq_set_si(four.val, 4, 1);
    rat_mul(&ac, a, c);
    rat_mul(&four_ac, &four, &ac);
    rat_sub(&d, &b2, &four_ac);

    if (rat_is_zero(&d)) {
        Rational neg_b, two_a, root;
        rat_init(&neg_b); rat_init(&two_a); rat_init(&root);
        rat_neg(&neg_b, b);
        Rational two = rat_from_int(2);
        rat_mul(&two_a, &two, a);
        rat_div(&root, &neg_b, &two_a);
        sol->rational_roots = realloc(sol->rational_roots,
            (sol->num_rational + 1) * sizeof(Rational));
        rat_init_set(&sol->rational_roots[sol->num_rational++], &root);
        rat_clear(&neg_b); rat_clear(&two_a); rat_clear(&root); rat_clear(&two);
    } else if (rat_is_negative(&d)) {
        Polynomial* irr = create_polynomial();
        add_term(irr, a, 2);
        add_term(irr, b, 1);
        add_term(irr, c, 0);
        sol->irreducible_factors = realloc(sol->irreducible_factors,
            (sol->num_irreducible + 1) * sizeof(Polynomial*));
        sol->irreducible_factors[sol->num_irreducible++] = irr;
    } else {
        Rational sqrt_d;
        rat_init(&sqrt_d);
        if (is_perfect_square_rat(&d, &sqrt_d)) {
            Rational neg_b, two_a, sum, diff, r1, r2;
            rat_init(&neg_b); rat_init(&two_a); rat_init(&sum); rat_init(&diff);
            rat_init(&r1); rat_init(&r2);
            rat_neg(&neg_b, b);
            Rational two = rat_from_int(2);
            rat_mul(&two_a, &two, a);
            rat_add(&sum, &neg_b, &sqrt_d);
            rat_div(&r1, &sum, &two_a);
            rat_sub(&diff, &neg_b, &sqrt_d);
            rat_div(&r2, &diff, &two_a);
            sol->rational_roots = realloc(sol->rational_roots,
                (sol->num_rational + 2) * sizeof(Rational));
            rat_init_set(&sol->rational_roots[sol->num_rational++], &r1);
            if (!rat_eq(&r1, &r2))
                rat_init_set(&sol->rational_roots[sol->num_rational++], &r2);
            rat_clear(&neg_b); rat_clear(&two_a); rat_clear(&sum); rat_clear(&diff);
            rat_clear(&r1); rat_clear(&r2); rat_clear(&two);
        } else {
            Rational neg_b, two_a;
            rat_init(&neg_b); rat_init(&two_a);
            rat_neg(&neg_b, b);
            Rational two = rat_from_int(2);
            rat_mul(&two_a, &two, a);
            double dd_val = sqrt(rat_to_double(&d));
            double bd = rat_to_double(&neg_b);
            double ad = rat_to_double(&two_a);
            sol->irrational_roots = realloc(sol->irrational_roots,
                (sol->num_irrational + 2) * sizeof(double));
            sol->irrational_roots[sol->num_irrational++] = (bd + dd_val) / ad;
            sol->irrational_roots[sol->num_irrational++] = (bd - dd_val) / ad;
            rat_clear(&neg_b); rat_clear(&two_a); rat_clear(&two);
        }
        rat_clear(&sqrt_d);
    }

    rat_clear(&b2); rat_clear(&four); rat_clear(&ac); rat_clear(&four_ac); rat_clear(&d);
}

static void get_coeff(Rational* result, const Polynomial* poly, int degree) {
    Term* t = poly->head;
    while (t) {
        if (t->exponent == degree) { rat_set(result, &t->coeff); return; }
        t = t->next;
    }
    mpq_set_si(result->val, 0, 1);
}

Solutions* solve_polynomial(const Polynomial* poly) {
    Solutions* sol = malloc(sizeof(Solutions));
    sol->rational_roots = NULL;
    sol->num_rational = 0;
    sol->irrational_roots = NULL;
    sol->num_irrational = 0;
    sol->irreducible_factors = NULL;
    sol->num_irreducible = 0;

    if (poly->head == NULL) return sol;

    Factorization* f = factorize(poly);
    if (!f) return sol;

    for (int i = 0; i < f->count; i++) {
        Polynomial* factor = f->factors[i];
        int deg = polynomial_degree(factor);
        int mult = f->multiplicities[i];

        if (deg == 1) {
            Rational a, b, root;
            rat_init(&a); rat_init(&b); rat_init(&root);
            get_coeff(&a, factor, 1);
            get_coeff(&b, factor, 0);
            Rational neg_b;
            rat_init(&neg_b);
            rat_neg(&neg_b, &b);
            rat_div(&root, &neg_b, &a);
            for (int m = 0; m < mult; m++) {
                sol->rational_roots = realloc(sol->rational_roots,
                    (sol->num_rational + 1) * sizeof(Rational));
                rat_init_set(&sol->rational_roots[sol->num_rational++], &root);
            }
            rat_clear(&a); rat_clear(&b); rat_clear(&root); rat_clear(&neg_b);
        } else if (deg == 2) {
            Rational a, b, c;
            rat_init(&a); rat_init(&b); rat_init(&c);
            get_coeff(&a, factor, 2);
            get_coeff(&b, factor, 1);
            get_coeff(&c, factor, 0);
            for (int m = 0; m < mult; m++)
                solve_quadratic(&a, &b, &c, sol);
            rat_clear(&a); rat_clear(&b); rat_clear(&c);
        } else {
            for (int m = 0; m < mult; m++) {
                sol->irreducible_factors = realloc(sol->irreducible_factors,
                    (sol->num_irreducible + 1) * sizeof(Polynomial*));
                sol->irreducible_factors[sol->num_irreducible++] = poly_copy(factor);
            }
        }
    }

    free_factorization(f);
    return sol;
}

char* solutions_to_string(const Solutions* sol) {
    size_t cap = 128, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';

    if (sol->num_rational == 0 && sol->num_irrational == 0 && sol->num_irreducible == 0) {
        strcpy(buf, "No solutions found");
        return buf;
    }

    for (int i = 0; i < sol->num_rational; i++) {
        char* rs = rat_to_string(&sol->rational_roots[i]);
        char line[128];
        snprintf(line, sizeof(line), "x = %s", rs);
        free(rs);
        size_t llen = strlen(line);
        while (len + llen + 3 > cap) { cap *= 2; buf = realloc(buf, cap); }
        if (len > 0) { strcat(buf, ", "); len += 2; }
        strcat(buf, line);
        len += llen;
    }

    for (int i = 0; i < sol->num_irrational; i++) {
        char line[128];
        snprintf(line, sizeof(line), "x ~ %.6f", sol->irrational_roots[i]);
        size_t llen = strlen(line);
        while (len + llen + 3 > cap) { cap *= 2; buf = realloc(buf, cap); }
        if (len > 0) { strcat(buf, ", "); len += 2; }
        strcat(buf, line);
        len += llen;
    }

    for (int i = 0; i < sol->num_irreducible; i++) {
        char* ps = poly_to_string(sol->irreducible_factors[i]);
        char line[256];
        snprintf(line, sizeof(line), "(%s) = 0 (no rational roots)", ps);
        free(ps);
        size_t llen = strlen(line);
        while (len + llen + 3 > cap) { cap *= 2; buf = realloc(buf, cap); }
        if (len > 0) { strcat(buf, ", "); len += 2; }
        strcat(buf, line);
        len += llen;
    }

    return buf;
}

void free_solutions(Solutions* sol) {
    if (!sol) return;
    for (int i = 0; i < sol->num_rational; i++)
        rat_clear(&sol->rational_roots[i]);
    free(sol->rational_roots);
    free(sol->irrational_roots);
    for (int i = 0; i < sol->num_irreducible; i++)
        free_polynomial(sol->irreducible_factors[i]);
    free(sol->irreducible_factors);
    free(sol);
}
