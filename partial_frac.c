#include "partial_frac.h"
#include "factor.h"
#include "matrix.h"
#include "strbuf.h"
#include <stdlib.h>

// content · Π_{j≠i} P_j^(m_j): the denominator with factor i removed entirely
static Polynomial* cofactor(const Factorization* F, int i) {
    Polynomial* c = poly_constant(&F->content);
    for (int j = 0; j < F->count; j++) {
        if (j == i) continue;
        Polynomial* pw = poly_pow(F->factors[j], (unsigned)F->multiplicities[j]);
        Polynomial* t = poly_multiply(c, pw);
        free_polynomial(c); free_polynomial(pw);
        c = t;
    }
    return c;
}

typedef struct {
    Polynomial** items;
    int count;
} PolyList;

// basis polynomials x^j · cofactor_i · P_i^(m-k) for k = 1..m, j = 0..deg P_i - 1,
// one per unknown coefficient, appended to the list
static void append_basis(const Factorization* F, int i, PolyList* basis) {
    Polynomial* co = cofactor(F, i);
    int m = F->multiplicities[i], d = polynomial_degree(F->factors[i]);
    Rational one = rat_one_val();
    for (int k = 1; k <= m; k++) {
        Polynomial* pw = poly_pow(F->factors[i], (unsigned)(m - k));
        Polynomial* base = poly_multiply(co, pw);
        for (int j = 0; j < d; j++) {
            Polynomial* xj = create_polynomial();
            add_term(xj, &one, j);
            basis->items[basis->count++] = poly_multiply(xj, base);
            free_polynomial(xj);
        }
        free_polynomial(pw); free_polynomial(base);
    }
    rat_clear(&one);
    free_polynomial(co);
}

// solve Σ u_k · basis_k = rem for the unknowns u (square system: deg den unknowns)
static LinSysResult* solve_coefficients(Polynomial** basis, int n, const Polynomial* rem) {
    Matrix* A = matrix_create(n, n);
    Vector* b = vector_create(n);
    for (int row = 0; row < n; row++) { // one equation per power of x
        int power = n - 1 - row;
        poly_coeff(&b->data[row], rem, power);
        for (int col = 0; col < n; col++) poly_coeff(&MAT_AT(A, row, col), basis[col], power);
    }
    LinSysResult* sol = linsys_solve(A, b);
    matrix_free(A);
    vector_free(b);
    return sol;
}

// read the solved unknowns back into per-power numerators
static void fill_term(PFTerm* t, const Vector* x, int* idx) {
    int d = polynomial_degree(t->factor);
    for (int k = 0; k < t->power; k++) {
        t->numerators[k] = create_polynomial();
        for (int j = 0; j < d; j++) add_term(t->numerators[k], &x->data[(*idx)++], j);
        if (t->coeffs) poly_coeff(&t->coeffs[k], t->numerators[k], 0);
    }
}

static void init_term(PFTerm* t, const Factorization* F, int i) {
    t->factor = poly_copy(F->factors[i]);
    t->power = F->multiplicities[i];
    t->numerators = calloc((size_t)t->power, sizeof(Polynomial*));
    t->coeffs = NULL;
    if (polynomial_degree(t->factor) == 1) {
        t->coeffs = malloc((size_t)t->power * sizeof(Rational));
        for (int k = 0; k < t->power; k++) rat_init(&t->coeffs[k]);
    }
}

PartialFractionDecomp* partial_fractions(const Polynomial* num, const Polynomial* den) {
    PartialFractionDecomp* pf = calloc(1, sizeof(PartialFractionDecomp));
    PolyDivResult dr = poly_divmod(num, den);
    if (!poly_is_zero(dr.quotient)) pf->poly_part = dr.quotient;
    else free_polynomial(dr.quotient);
    if (poly_is_zero(dr.remainder) || polynomial_degree(den) <= 0) {
        free_polynomial(dr.remainder);
        return pf;
    }
    Factorization* F = factorize(den);
    int n = polynomial_degree(den);
    PolyList basis = { malloc((size_t)n * sizeof(Polynomial*)), 0 };
    for (int i = 0; i < F->count; i++) append_basis(F, i, &basis);
    LinSysResult* sol = solve_coefficients(basis.items, n, dr.remainder);
    pf->num_terms = F->count;
    pf->terms = calloc((size_t)F->count, sizeof(PFTerm));
    int idx = 0;
    for (int i = 0; i < F->count; i++) {
        init_term(&pf->terms[i], F, i);
        fill_term(&pf->terms[i], sol->solution, &idx); // the system is always uniquely solvable
    }
    for (int i = 0; i < n; i++) free_polynomial(basis.items[i]);
    free(basis.items);
    linsys_result_free(sol);
    free_factorization(F);
    free_polynomial(dr.remainder);
    return pf;
}

void free_pf_decomp(PartialFractionDecomp* pf) {
    if (!pf) return;
    free_polynomial(pf->poly_part);
    for (int i = 0; i < pf->num_terms; i++) {
        PFTerm* t = &pf->terms[i];
        for (int k = 0; k < t->power; k++) {
            free_polynomial(t->numerators[k]);
            if (t->coeffs) rat_clear(&t->coeffs[k]);
        }
        free(t->numerators);
        free(t->coeffs);
        free_polynomial(t->factor);
    }
    free(pf->terms);
    free(pf);
}

typedef struct {
    const Polynomial* N; // numerator
    const char* fs;      // printed factor P
    int k;               // power of P
} Fraction;

// one summand N/(P)^k, with constant numerators folded into the sign
static void append_fraction(StrBuf* sb, const Fraction* fr) {
    const Polynomial* N = fr->N;
    bool first = sb->len == 0;
    if (polynomial_degree(N) == 0) {
        bool neg = rat_is_negative(&N->head->coeff);
        Rational a;
        rat_init(&a);
        rat_abs(&a, &N->head->coeff);
        sb_append(sb, first ? (neg ? "-" : "") : (neg ? " - " : " + "));
        sb_append_owned(sb, rat_to_string(&a));
        rat_clear(&a);
    } else {
        sb_append(sb, first ? "(" : " + (");
        sb_append_owned(sb, poly_to_string(N));
        sb_append(sb, ")");
    }
    sb_appendf(sb, "/(%s)", fr->fs);
    if (fr->k > 1) sb_appendf(sb, "^%d", fr->k);
}

char* pf_to_string(const PartialFractionDecomp* pf) {
    StrBuf sb = sb_new();
    if (pf->poly_part) sb_append_owned(&sb, poly_to_string(pf->poly_part));
    for (int i = 0; i < pf->num_terms; i++) {
        char* fs = poly_to_string(pf->terms[i].factor);
        for (int k = 0; k < pf->terms[i].power; k++)
            if (!poly_is_zero(pf->terms[i].numerators[k])) {
                Fraction fr = { pf->terms[i].numerators[k], fs, k + 1 };
                append_fraction(&sb, &fr);
            }
        free(fs);
    }
    if (sb.len == 0) sb_append(&sb, "0");
    return sb_finish(&sb);
}
