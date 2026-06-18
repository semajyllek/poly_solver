#include "solver.h"
#include "factor.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// check if a rational is a perfect square, return its square root
static bool is_perfect_square_rat(Rational r, Rational* root) {
    if (rat_is_negative(r)) return false;
    if (rat_is_zero(r)) { *root = rat_zero(); return true; }

    double nd = sqrt((double)r.num);
    double dd = sqrt((double)r.den);
    rat_int_t ni = (rat_int_t)round(nd);
    rat_int_t di = (rat_int_t)round(dd);
    if (ni * ni == r.num && di * di == r.den) {
        *root = rat_create(ni, di);
        return true;
    }
    return false;
}

// solve ax^2 + bx + c = 0
static void solve_quadratic(Rational a, Rational b, Rational c, Solutions* sol) {
    // discriminant D = b^2 - 4ac
    Rational d = rat_sub(rat_mul(b, b), rat_mul(rat_from_int(4), rat_mul(a, c)));

    if (rat_is_zero(d)) {
        // one repeated root: -b / (2a)
        Rational root = rat_div(rat_neg(b), rat_mul(rat_from_int(2), a));
        sol->rational_roots = realloc(sol->rational_roots,
            (sol->num_rational + 1) * sizeof(Rational));
        sol->rational_roots[sol->num_rational++] = root;
        return;
    }

    if (rat_is_negative(d)) {
        // no real roots — record the irreducible quadratic
        Polynomial* irr = create_polynomial();
        add_term(irr, a, 2);
        add_term(irr, b, 1);
        add_term(irr, c, 0);
        sol->irreducible_factors = realloc(sol->irreducible_factors,
            (sol->num_irreducible + 1) * sizeof(Polynomial*));
        sol->irreducible_factors[sol->num_irreducible++] = irr;
        return;
    }

    // D > 0
    Rational sqrt_d;
    if (is_perfect_square_rat(d, &sqrt_d)) {
        // exact rational roots
        Rational two_a = rat_mul(rat_from_int(2), a);
        Rational r1 = rat_div(rat_add(rat_neg(b), sqrt_d), two_a);
        Rational r2 = rat_div(rat_sub(rat_neg(b), sqrt_d), two_a);
        sol->rational_roots = realloc(sol->rational_roots,
            (sol->num_rational + 2) * sizeof(Rational));
        sol->rational_roots[sol->num_rational++] = r1;
        if (!rat_eq(r1, r2))
            sol->rational_roots[sol->num_rational++] = r2;
    } else {
        // irrational roots — return as doubles
        double dd = sqrt(rat_to_double(d));
        double bd = rat_to_double(rat_neg(b));
        double ad = rat_to_double(rat_mul(rat_from_int(2), a));
        sol->irrational_roots = realloc(sol->irrational_roots,
            (sol->num_irrational + 2) * sizeof(double));
        sol->irrational_roots[sol->num_irrational++] = (bd + dd) / ad;
        sol->irrational_roots[sol->num_irrational++] = (bd - dd) / ad;
    }
}

// get coefficient at given degree from a polynomial
static Rational get_coeff(const Polynomial* poly, int degree) {
    Term* t = poly->head;
    while (t) {
        if (t->exponent == degree) return t->coeff;
        t = t->next;
    }
    return rat_zero();
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

    // factor first, then solve each factor
    Factorization* f = factorize(poly);
    if (!f) return sol;

    for (int i = 0; i < f->count; i++) {
        Polynomial* factor = f->factors[i];
        int deg = polynomial_degree(factor);
        int mult = f->multiplicities[i];

        if (deg == 1) {
            Rational a = get_coeff(factor, 1);
            Rational b = get_coeff(factor, 0);
            // root: -b/a, with multiplicity
            Rational root = rat_neg(rat_div(b, a));
            for (int m = 0; m < mult; m++) {
                sol->rational_roots = realloc(sol->rational_roots,
                    (sol->num_rational + 1) * sizeof(Rational));
                sol->rational_roots[sol->num_rational++] = root;
            }
        } else if (deg == 2) {
            Rational a = get_coeff(factor, 2);
            Rational b = get_coeff(factor, 1);
            Rational c = get_coeff(factor, 0);
            for (int m = 0; m < mult; m++)
                solve_quadratic(a, b, c, sol);
        } else {
            // degree > 2 with no rational roots: report as irreducible
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
        char* rs = rat_to_string(sol->rational_roots[i]);
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

    if (sol->num_irreducible > 0) {
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
    }

    return buf;
}

void free_solutions(Solutions* sol) {
    if (!sol) return;
    free(sol->rational_roots);
    free(sol->irrational_roots);
    for (int i = 0; i < sol->num_irreducible; i++)
        free_polynomial(sol->irreducible_factors[i]);
    free(sol->irreducible_factors);
    free(sol);
}
