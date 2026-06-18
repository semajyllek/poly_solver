#include "partial_frac.h"
#include "factor.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// get coefficient of x^degree from a polynomial
static Rational get_coeff(const Polynomial* poly, int degree) {
    Term* t = poly->head;
    while (t) {
        if (t->exponent == degree) return t->coeff;
        t = t->next;
    }
    return rat_zero();
}

PartialFractionDecomp* partial_fractions(const Polynomial* num, const Polynomial* den) {
    PartialFractionDecomp* pf = malloc(sizeof(PartialFractionDecomp));
    pf->poly_part = NULL;
    pf->terms = NULL;
    pf->num_terms = 0;

    // step 1: ensure proper fraction via long division
    Polynomial* remainder;
    if (polynomial_degree(num) >= polynomial_degree(den)) {
        PolyDivResult dr = poly_divmod(num, den);
        pf->poly_part = dr.quotient;
        remainder = dr.remainder;
    } else {
        remainder = poly_copy(num);
    }

    // trivial case: remainder is zero
    if (poly_is_zero(remainder)) {
        free_polynomial(remainder);
        return pf;
    }

    // step 2: factor denominator
    Factorization* facts = factorize(den);
    if (!facts || facts->count == 0) {
        free_factorization(facts);
        free_polynomial(remainder);
        return pf;
    }

    // step 3: count unknowns and allocate terms
    int num_unknowns = 0;
    pf->num_terms = facts->count;
    pf->terms = calloc(facts->count, sizeof(PFTerm));

    for (int i = 0; i < facts->count; i++) {
        int deg = polynomial_degree(facts->factors[i]);
        int mult = facts->multiplicities[i];
        pf->terms[i].factor = poly_copy(facts->factors[i]);
        pf->terms[i].power = mult;

        if (deg == 1) {
            pf->terms[i].is_quadratic = false;
            pf->terms[i].coeffs = calloc(mult, sizeof(Rational));
            pf->terms[i].quad_coeffs = NULL;
            num_unknowns += mult;
        } else {
            // irreducible quadratic or higher — store as quadratic
            pf->terms[i].is_quadratic = true;
            pf->terms[i].coeffs = NULL;
            pf->terms[i].quad_coeffs = calloc(2 * mult, sizeof(Rational));
            num_unknowns += 2 * mult;
        }
    }

    // step 4: build and solve the linear system
    // the system has num_unknowns columns and deg(den) rows (after content extraction)
    int den_degree = polynomial_degree(den);
    // the content was factored out by factorize(), so we work with the original den

    // for each unknown, compute the polynomial it multiplies when clearing denominators
    // unknown_i corresponds to A_i/(factor^k)
    // clearing denominators: A_i * (den / factor^k)
    Polynomial** basis_polys = malloc(num_unknowns * sizeof(Polynomial*));
    int unknown_idx = 0;

    for (int i = 0; i < facts->count; i++) {
        int mult = facts->multiplicities[i];
        int deg = polynomial_degree(facts->factors[i]);

        // compute den / factor^mult (the cofactor for this factor group)
        Polynomial* cofactor = create_polynomial();
        add_term(cofactor, facts->content, 0);
        for (int j = 0; j < facts->count; j++) {
            if (j == i) continue;
            for (int m = 0; m < facts->multiplicities[j]; m++) {
                Polynomial* tmp = poly_multiply(cofactor, facts->factors[j]);
                free_polynomial(cofactor);
                cofactor = tmp;
            }
        }

        if (deg == 1) {
            // for each power k = 1..mult: unknown * factor^(mult-k) * cofactor
            for (int k = 1; k <= mult; k++) {
                // multiply cofactor by factor^(mult-k)
                Polynomial* bp = poly_copy(cofactor);
                for (int p = 0; p < mult - k; p++) {
                    Polynomial* tmp = poly_multiply(bp, facts->factors[i]);
                    free_polynomial(bp);
                    bp = tmp;
                }
                basis_polys[unknown_idx++] = bp;
            }
        } else {
            // quadratic: for each power k = 1..mult, two unknowns (Bx + C)
            for (int k = 1; k <= mult; k++) {
                Polynomial* bp = poly_copy(cofactor);
                for (int p = 0; p < mult - k; p++) {
                    Polynomial* tmp = poly_multiply(bp, facts->factors[i]);
                    free_polynomial(bp);
                    bp = tmp;
                }
                // B * x * bp
                Polynomial* x_poly = create_polynomial();
                add_term(x_poly, rat_one(), 1);
                Polynomial* x_bp = poly_multiply(x_poly, bp);
                basis_polys[unknown_idx++] = x_bp;
                // C * bp
                basis_polys[unknown_idx++] = poly_copy(bp);
                free_polynomial(x_poly);
                free_polynomial(bp);
            }
        }
        free_polynomial(cofactor);
    }

    // build the matrix: rows = den_degree, cols = num_unknowns
    Matrix* A = matrix_create(den_degree, num_unknowns);
    Vector* b = vector_create(den_degree);

    for (int row = 0; row < den_degree; row++) {
        int power = den_degree - 1 - row;
        b->data[row] = get_coeff(remainder, power);
        for (int col = 0; col < num_unknowns; col++)
            MAT_AT(A, row, col) = get_coeff(basis_polys[col], power);
    }

    LinSysResult* sol = linsys_solve(A, b);

    if (sol->status == LINSYS_UNIQUE) {
        // distribute solution back to PFTerms
        unknown_idx = 0;
        for (int i = 0; i < pf->num_terms; i++) {
            int mult = pf->terms[i].power;
            if (!pf->terms[i].is_quadratic) {
                for (int k = 0; k < mult; k++)
                    pf->terms[i].coeffs[k] = sol->solution->data[unknown_idx++];
            } else {
                for (int k = 0; k < mult; k++) {
                    pf->terms[i].quad_coeffs[2 * k] = sol->solution->data[unknown_idx++];
                    pf->terms[i].quad_coeffs[2 * k + 1] = sol->solution->data[unknown_idx++];
                }
            }
        }
    }

    // cleanup
    for (int i = 0; i < num_unknowns; i++)
        free_polynomial(basis_polys[i]);
    free(basis_polys);
    linsys_result_free(sol);
    matrix_free(A);
    vector_free(b);
    free_polynomial(remainder);
    free_factorization(facts);

    return pf;
}

void free_pf_decomp(PartialFractionDecomp* pf) {
    if (!pf) return;
    free_polynomial(pf->poly_part);
    for (int i = 0; i < pf->num_terms; i++) {
        free(pf->terms[i].coeffs);
        free(pf->terms[i].quad_coeffs);
        free_polynomial(pf->terms[i].factor);
    }
    free(pf->terms);
    free(pf);
}

char* pf_to_string(const PartialFractionDecomp* pf) {
    size_t cap = 256, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';

    bool has_content = false;

    // polynomial part
    if (pf->poly_part && !poly_is_zero(pf->poly_part)) {
        char* ps = poly_to_string(pf->poly_part);
        while (len + strlen(ps) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, ps);
        len += strlen(ps);
        free(ps);
        has_content = true;
    }

    for (int i = 0; i < pf->num_terms; i++) {
        PFTerm* t = &pf->terms[i];
        char* fs = poly_to_string(t->factor);

        if (!t->is_quadratic) {
            for (int k = 0; k < t->power; k++) {
                Rational coeff = t->coeffs[k];
                if (rat_is_zero(coeff)) continue;

                char term_buf[256];
                char* cs = rat_to_string(rat_abs(coeff));
                int pw = k + 1; // power of the factor in denominator

                if (pw == 1)
                    snprintf(term_buf, sizeof(term_buf), "%s/(%s)", cs, fs);
                else
                    snprintf(term_buf, sizeof(term_buf), "%s/(%s)^%d", cs, fs, pw);
                free(cs);

                // separator
                if (has_content) {
                    if (rat_is_negative(coeff)) {
                        while (len + 4 > cap) { cap *= 2; buf = realloc(buf, cap); }
                        strcat(buf, " - "); len += 3;
                    } else {
                        while (len + 4 > cap) { cap *= 2; buf = realloc(buf, cap); }
                        strcat(buf, " + "); len += 3;
                    }
                } else if (rat_is_negative(coeff)) {
                    while (len + 2 > cap) { cap *= 2; buf = realloc(buf, cap); }
                    strcat(buf, "-"); len += 1;
                }

                size_t tlen = strlen(term_buf);
                while (len + tlen + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
                strcat(buf, term_buf);
                len += tlen;
                has_content = true;
            }
        } else {
            for (int k = 0; k < t->power; k++) {
                Rational B = t->quad_coeffs[2 * k];
                Rational C = t->quad_coeffs[2 * k + 1];
                if (rat_is_zero(B) && rat_is_zero(C)) continue;

                char term_buf[256];
                char* bs = rat_to_string(B);
                char* ccs = rat_to_string(C);
                int pw = k + 1;
                if (pw == 1)
                    snprintf(term_buf, sizeof(term_buf), "(%sx + %s)/(%s)", bs, ccs, fs);
                else
                    snprintf(term_buf, sizeof(term_buf), "(%sx + %s)/(%s)^%d", bs, ccs, fs, pw);
                free(bs);
                free(ccs);

                if (has_content) {
                    while (len + 4 > cap) { cap *= 2; buf = realloc(buf, cap); }
                    strcat(buf, " + "); len += 3;
                }

                size_t tlen = strlen(term_buf);
                while (len + tlen + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
                strcat(buf, term_buf);
                len += tlen;
                has_content = true;
            }
        }
        free(fs);
    }

    if (!has_content) { strcpy(buf, "0"); }
    return buf;
}
