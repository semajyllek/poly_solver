#include "partial_frac.h"
#include "factor.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void pf_get_coeff(Rational* result, const Polynomial* poly, int degree) {
    Term* t = poly->head;
    while (t) {
        if (t->exponent == degree) { rat_set(result, &t->coeff); return; }
        t = t->next;
    }
    mpq_set_si(result->val, 0, 1);
}

PartialFractionDecomp* partial_fractions(const Polynomial* num, const Polynomial* den) {
    PartialFractionDecomp* pf = malloc(sizeof(PartialFractionDecomp));
    pf->poly_part = NULL;
    pf->terms = NULL;
    pf->num_terms = 0;

    Polynomial* remainder;
    if (polynomial_degree(num) >= polynomial_degree(den)) {
        PolyDivResult dr = poly_divmod(num, den);
        pf->poly_part = dr.quotient;
        remainder = dr.remainder;
    } else {
        remainder = poly_copy(num);
    }

    if (poly_is_zero(remainder)) {
        free_polynomial(remainder);
        return pf;
    }

    Factorization* facts = factorize(den);
    if (!facts || facts->count == 0) {
        free_factorization(facts);
        free_polynomial(remainder);
        return pf;
    }

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
            pf->terms[i].coeffs = malloc(mult * sizeof(Rational));
            for (int k = 0; k < mult; k++) rat_init(&pf->terms[i].coeffs[k]);
            pf->terms[i].quad_coeffs = NULL;
            num_unknowns += mult;
        } else {
            pf->terms[i].is_quadratic = true;
            pf->terms[i].coeffs = NULL;
            pf->terms[i].quad_coeffs = malloc(2 * mult * sizeof(Rational));
            for (int k = 0; k < 2 * mult; k++) rat_init(&pf->terms[i].quad_coeffs[k]);
            num_unknowns += 2 * mult;
        }
    }

    int den_degree = polynomial_degree(den);

    Polynomial** basis_polys = malloc(num_unknowns * sizeof(Polynomial*));
    int unknown_idx = 0;

    for (int i = 0; i < facts->count; i++) {
        int mult = facts->multiplicities[i];
        int deg = polynomial_degree(facts->factors[i]);

        Polynomial* cofactor = create_polynomial();
        add_term(cofactor, &facts->content, 0);
        for (int j = 0; j < facts->count; j++) {
            if (j == i) continue;
            for (int m = 0; m < facts->multiplicities[j]; m++) {
                Polynomial* tmp = poly_multiply(cofactor, facts->factors[j]);
                free_polynomial(cofactor);
                cofactor = tmp;
            }
        }

        if (deg == 1) {
            for (int k = 1; k <= mult; k++) {
                Polynomial* bp = poly_copy(cofactor);
                for (int p = 0; p < mult - k; p++) {
                    Polynomial* tmp = poly_multiply(bp, facts->factors[i]);
                    free_polynomial(bp);
                    bp = tmp;
                }
                basis_polys[unknown_idx++] = bp;
            }
        } else {
            for (int k = 1; k <= mult; k++) {
                Polynomial* bp = poly_copy(cofactor);
                for (int p = 0; p < mult - k; p++) {
                    Polynomial* tmp = poly_multiply(bp, facts->factors[i]);
                    free_polynomial(bp);
                    bp = tmp;
                }
                Polynomial* x_poly = create_polynomial();
                Rational one_r = rat_one_val();
                add_term(x_poly, &one_r, 1);
                rat_clear(&one_r);
                Polynomial* x_bp = poly_multiply(x_poly, bp);
                basis_polys[unknown_idx++] = x_bp;
                basis_polys[unknown_idx++] = poly_copy(bp);
                free_polynomial(x_poly);
                free_polynomial(bp);
            }
        }
        free_polynomial(cofactor);
    }

    Matrix* A = matrix_create(den_degree, num_unknowns);
    Vector* b = vector_create(den_degree);

    Rational tmp;
    rat_init(&tmp);
    for (int row = 0; row < den_degree; row++) {
        int power = den_degree - 1 - row;
        pf_get_coeff(&b->data[row], remainder, power);
        for (int col = 0; col < num_unknowns; col++) {
            pf_get_coeff(&tmp, basis_polys[col], power);
            rat_set(&MAT_AT(A, row, col), &tmp);
        }
    }
    rat_clear(&tmp);

    LinSysResult* sol = linsys_solve(A, b);

    if (sol->status == LINSYS_UNIQUE) {
        unknown_idx = 0;
        for (int i = 0; i < pf->num_terms; i++) {
            int mult = pf->terms[i].power;
            if (!pf->terms[i].is_quadratic) {
                for (int k = 0; k < mult; k++)
                    rat_set(&pf->terms[i].coeffs[k], &sol->solution->data[unknown_idx++]);
            } else {
                for (int k = 0; k < mult; k++) {
                    rat_set(&pf->terms[i].quad_coeffs[2*k], &sol->solution->data[unknown_idx++]);
                    rat_set(&pf->terms[i].quad_coeffs[2*k+1], &sol->solution->data[unknown_idx++]);
                }
            }
        }
    }

    for (int i = 0; i < num_unknowns; i++) free_polynomial(basis_polys[i]);
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
        if (pf->terms[i].coeffs) {
            for (int k = 0; k < pf->terms[i].power; k++)
                rat_clear(&pf->terms[i].coeffs[k]);
            free(pf->terms[i].coeffs);
        }
        if (pf->terms[i].quad_coeffs) {
            for (int k = 0; k < 2 * pf->terms[i].power; k++)
                rat_clear(&pf->terms[i].quad_coeffs[k]);
            free(pf->terms[i].quad_coeffs);
        }
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

    if (pf->poly_part && !poly_is_zero(pf->poly_part)) {
        char* ps = poly_to_string(pf->poly_part);
        while (len + strlen(ps) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, ps); len += strlen(ps);
        free(ps);
        has_content = true;
    }

    for (int i = 0; i < pf->num_terms; i++) {
        PFTerm* t = &pf->terms[i];
        char* fs = poly_to_string(t->factor);

        if (!t->is_quadratic) {
            for (int k = 0; k < t->power; k++) {
                if (rat_is_zero(&t->coeffs[k])) continue;

                Rational abs_c;
                rat_init(&abs_c);
                rat_abs(&abs_c, &t->coeffs[k]);
                char* cs = rat_to_string(&abs_c);
                rat_clear(&abs_c);
                int pw = k + 1;
                char term_buf[256];
                if (pw == 1) snprintf(term_buf, sizeof(term_buf), "%s/(%s)", cs, fs);
                else snprintf(term_buf, sizeof(term_buf), "%s/(%s)^%d", cs, fs, pw);
                free(cs);

                if (has_content) {
                    const char* sep = rat_is_negative(&t->coeffs[k]) ? " - " : " + ";
                    while (len + 4 > cap) { cap *= 2; buf = realloc(buf, cap); }
                    strcat(buf, sep); len += 3;
                } else if (rat_is_negative(&t->coeffs[k])) {
                    while (len + 2 > cap) { cap *= 2; buf = realloc(buf, cap); }
                    strcat(buf, "-"); len += 1;
                }

                size_t tlen = strlen(term_buf);
                while (len + tlen + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
                strcat(buf, term_buf); len += tlen;
                has_content = true;
            }
        } else {
            for (int k = 0; k < t->power; k++) {
                if (rat_is_zero(&t->quad_coeffs[2*k]) && rat_is_zero(&t->quad_coeffs[2*k+1])) continue;
                char* bs = rat_to_string(&t->quad_coeffs[2*k]);
                char* ccs = rat_to_string(&t->quad_coeffs[2*k+1]);
                int pw = k + 1;
                char term_buf[256];
                if (pw == 1) snprintf(term_buf, sizeof(term_buf), "(%sx + %s)/(%s)", bs, ccs, fs);
                else snprintf(term_buf, sizeof(term_buf), "(%sx + %s)/(%s)^%d", bs, ccs, fs, pw);
                free(bs); free(ccs);
                if (has_content) {
                    while (len + 4 > cap) { cap *= 2; buf = realloc(buf, cap); }
                    strcat(buf, " + "); len += 3;
                }
                size_t tlen = strlen(term_buf);
                while (len + tlen + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
                strcat(buf, term_buf); len += tlen;
                has_content = true;
            }
        }
        free(fs);
    }

    if (!has_content) strcpy(buf, "0");
    return buf;
}
