#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

Matrix* matrix_create(int rows, int cols) {
    Matrix* m = malloc(sizeof(Matrix));
    m->rows = rows;
    m->cols = cols;
    m->data = malloc(rows * cols * sizeof(Rational));
    for (int i = 0; i < rows * cols; i++)
        rat_init(&m->data[i]);
    return m;
}

Matrix* matrix_copy(const Matrix* m) {
    Matrix* copy = malloc(sizeof(Matrix));
    copy->rows = m->rows;
    copy->cols = m->cols;
    int n = m->rows * m->cols;
    copy->data = malloc(n * sizeof(Rational));
    for (int i = 0; i < n; i++)
        rat_init_set(&copy->data[i], &m->data[i]);
    return copy;
}

void matrix_free(Matrix* m) {
    if (!m) return;
    for (int i = 0; i < m->rows * m->cols; i++)
        rat_clear(&m->data[i]);
    free(m->data);
    free(m);
}

Vector* vector_create(int n) {
    Vector* v = malloc(sizeof(Vector));
    v->n = n;
    v->data = malloc(n * sizeof(Rational));
    for (int i = 0; i < n; i++)
        rat_init(&v->data[i]);
    return v;
}

Vector* vector_copy(const Vector* v) {
    Vector* copy = malloc(sizeof(Vector));
    copy->n = v->n;
    copy->data = malloc(v->n * sizeof(Rational));
    for (int i = 0; i < v->n; i++)
        rat_init_set(&copy->data[i], &v->data[i]);
    return copy;
}

void vector_free(Vector* v) {
    if (!v) return;
    for (int i = 0; i < v->n; i++)
        rat_clear(&v->data[i]);
    free(v->data);
    free(v);
}

static int gauss_eliminate(Matrix* m, int augmented_cols) {
    int pivot_row = 0, pivot_col = 0, swaps = 0;
    int coeff_cols = m->cols - augmented_cols;

    Rational tmp, factor, prod;
    rat_init(&tmp); rat_init(&factor); rat_init(&prod);

    while (pivot_row < m->rows && pivot_col < coeff_cols) {
        int found = -1;
        for (int r = pivot_row; r < m->rows; r++) {
            if (!rat_is_zero(&MAT_AT(m, r, pivot_col))) { found = r; break; }
        }
        if (found == -1) { pivot_col++; continue; }

        if (found != pivot_row) {
            for (int c = 0; c < m->cols; c++) {
                rat_set(&tmp, &MAT_AT(m, pivot_row, c));
                rat_set(&MAT_AT(m, pivot_row, c), &MAT_AT(m, found, c));
                rat_set(&MAT_AT(m, found, c), &tmp);
            }
            swaps++;
        }

        // scale pivot row
        Rational pivot_val;
        rat_init_set(&pivot_val, &MAT_AT(m, pivot_row, pivot_col));
        for (int c = pivot_col; c < m->cols; c++) {
            rat_div(&tmp, &MAT_AT(m, pivot_row, c), &pivot_val);
            rat_set(&MAT_AT(m, pivot_row, c), &tmp);
        }
        rat_clear(&pivot_val);

        // eliminate
        for (int r = 0; r < m->rows; r++) {
            if (r == pivot_row) continue;
            if (rat_is_zero(&MAT_AT(m, r, pivot_col))) continue;
            rat_set(&factor, &MAT_AT(m, r, pivot_col));
            for (int c = pivot_col; c < m->cols; c++) {
                rat_mul(&prod, &factor, &MAT_AT(m, pivot_row, c));
                rat_sub(&MAT_AT(m, r, c), &MAT_AT(m, r, c), &prod);
            }
        }

        pivot_row++; pivot_col++;
    }

    rat_clear(&tmp); rat_clear(&factor); rat_clear(&prod);
    return swaps;
}

Matrix* matrix_rref(const Matrix* m) {
    Matrix* result = matrix_copy(m);
    gauss_eliminate(result, 0);
    return result;
}

int matrix_rank(const Matrix* m) {
    Matrix* rref = matrix_rref(m);
    int rank = 0;
    for (int r = 0; r < rref->rows; r++) {
        for (int c = 0; c < rref->cols; c++) {
            if (!rat_is_zero(&MAT_AT(rref, r, c))) { rank++; break; }
        }
    }
    matrix_free(rref);
    return rank;
}

void matrix_det(Rational* result, const Matrix* m) {
    assert(m->rows == m->cols);
    int n = m->rows;
    Matrix* work = matrix_copy(m);
    mpq_set_si(result->val, 1, 1);
    int swaps = 0;

    Rational tmp, factor_r, prod;
    rat_init(&tmp); rat_init(&factor_r); rat_init(&prod);

    for (int col = 0; col < n; col++) {
        int found = -1;
        for (int r = col; r < n; r++) {
            if (!rat_is_zero(&MAT_AT(work, r, col))) { found = r; break; }
        }
        if (found == -1) {
            mpq_set_si(result->val, 0, 1);
            matrix_free(work);
            rat_clear(&tmp); rat_clear(&factor_r); rat_clear(&prod);
            return;
        }
        if (found != col) {
            for (int c = 0; c < n; c++) {
                rat_set(&tmp, &MAT_AT(work, col, c));
                rat_set(&MAT_AT(work, col, c), &MAT_AT(work, found, c));
                rat_set(&MAT_AT(work, found, c), &tmp);
            }
            swaps++;
        }

        rat_mul(result, result, &MAT_AT(work, col, col));

        for (int r = col + 1; r < n; r++) {
            rat_div(&factor_r, &MAT_AT(work, r, col), &MAT_AT(work, col, col));
            for (int c = col; c < n; c++) {
                rat_mul(&prod, &factor_r, &MAT_AT(work, col, c));
                rat_sub(&MAT_AT(work, r, c), &MAT_AT(work, r, c), &prod);
            }
        }
    }

    if (swaps % 2 == 1) rat_neg(result, result);
    matrix_free(work);
    rat_clear(&tmp); rat_clear(&factor_r); rat_clear(&prod);
}

Matrix* matrix_inverse(const Matrix* m) {
    assert(m->rows == m->cols);
    int n = m->rows;
    Matrix* aug = matrix_create(n, 2 * n);
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++)
            rat_set(&MAT_AT(aug, r, c), &MAT_AT(m, r, c));
        mpq_set_si(MAT_AT(aug, r, n + r).val, 1, 1);
    }
    gauss_eliminate(aug, n);

    Rational one_r = rat_one_val();
    Rational zero_r = rat_zero_val();
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            const Rational* expected = (r == c) ? &one_r : &zero_r;
            if (!rat_eq(&MAT_AT(aug, r, c), expected)) {
                matrix_free(aug);
                rat_clear(&one_r); rat_clear(&zero_r);
                return NULL;
            }
        }
    }
    rat_clear(&one_r); rat_clear(&zero_r);

    Matrix* inv = matrix_create(n, n);
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
            rat_set(&MAT_AT(inv, r, c), &MAT_AT(aug, r, n + c));

    matrix_free(aug);
    return inv;
}

LinSysResult* linsys_solve(const Matrix* A, const Vector* b) {
    assert(A->rows == b->n);
    LinSysResult* res = malloc(sizeof(LinSysResult));
    res->solution = NULL;
    res->rank = 0;

    Matrix* aug = matrix_create(A->rows, A->cols + 1);
    for (int r = 0; r < A->rows; r++) {
        for (int c = 0; c < A->cols; c++)
            rat_set(&MAT_AT(aug, r, c), &MAT_AT(A, r, c));
        rat_set(&MAT_AT(aug, r, A->cols), &b->data[r]);
    }
    gauss_eliminate(aug, 1);

    int rank = 0;
    for (int r = 0; r < aug->rows; r++) {
        for (int c = 0; c < A->cols; c++) {
            if (!rat_is_zero(&MAT_AT(aug, r, c))) { rank++; break; }
        }
    }
    res->rank = rank;

    for (int r = 0; r < aug->rows; r++) {
        bool all_zero = true;
        for (int c = 0; c < A->cols; c++) {
            if (!rat_is_zero(&MAT_AT(aug, r, c))) { all_zero = false; break; }
        }
        if (all_zero && !rat_is_zero(&MAT_AT(aug, r, A->cols))) {
            res->status = LINSYS_INCONSISTENT;
            matrix_free(aug);
            return res;
        }
    }

    if (rank < A->cols) {
        res->status = LINSYS_INFINITE;
        matrix_free(aug);
        return res;
    }

    res->status = LINSYS_UNIQUE;
    res->solution = vector_create(A->cols);
    for (int i = 0; i < A->cols; i++)
        rat_set(&res->solution->data[i], &MAT_AT(aug, i, A->cols));

    matrix_free(aug);
    return res;
}

// --- output ---

char* matrix_to_string(const Matrix* m) {
    size_t cap = 128, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';
    for (int r = 0; r < m->rows; r++) {
        const char* rs = (r == 0) ? "[[" : " [";
        while (len + strlen(rs) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, rs); len += strlen(rs);
        for (int c = 0; c < m->cols; c++) {
            char* vs = rat_to_string(&MAT_AT(m, r, c));
            const char* sep = (c < m->cols - 1) ? ", " : "";
            while (len + strlen(vs) + strlen(sep) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
            strcat(buf, vs); strcat(buf, sep);
            len += strlen(vs) + strlen(sep);
            free(vs);
        }
        const char* re = (r < m->rows - 1) ? "]\n" : "]]";
        while (len + strlen(re) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, re); len += strlen(re);
    }
    return buf;
}

char* vector_to_string(const Vector* v) {
    size_t cap = 64, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';
    strcat(buf, "["); len++;
    for (int i = 0; i < v->n; i++) {
        char* rs = rat_to_string(&v->data[i]);
        const char* sep = (i < v->n - 1) ? ", " : "";
        while (len + strlen(rs) + strlen(sep) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, rs); strcat(buf, sep);
        len += strlen(rs) + strlen(sep);
        free(rs);
    }
    while (len + 2 > cap) { cap *= 2; buf = realloc(buf, cap); }
    strcat(buf, "]");
    return buf;
}

char* linsys_result_to_string(const LinSysResult* r) {
    size_t cap = 128;
    char* buf = malloc(cap);
    if (r->status == LINSYS_INCONSISTENT) {
        snprintf(buf, cap, "No solution (inconsistent system, rank = %d)", r->rank);
        return buf;
    }
    if (r->status == LINSYS_INFINITE) {
        snprintf(buf, cap, "Infinitely many solutions (rank = %d, %d variables)", r->rank, r->rank);
        return buf;
    }
    size_t len = 0;
    buf[0] = '\0';
    for (int i = 0; i < r->solution->n; i++) {
        char* rs = rat_to_string(&r->solution->data[i]);
        char entry[128];
        snprintf(entry, sizeof(entry), "x%d = %s", i + 1, rs);
        free(rs);
        const char* sep = (i < r->solution->n - 1) ? ", " : "";
        size_t elen = strlen(entry);
        while (len + elen + strlen(sep) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, entry); strcat(buf, sep);
        len += elen + strlen(sep);
    }
    return buf;
}

void linsys_result_free(LinSysResult* r) {
    if (!r) return;
    vector_free(r->solution);
    free(r);
}

// --- CLI parsing ---

static void skip_ws(const char** s) { while (isspace(**s)) (*s)++; }

static bool parse_rational_lit(const char** s, Rational* out) {
    skip_ws(s);
    char* end;
    long num = strtol(*s, &end, 10);
    if (end == *s) return false;
    *s = end;
    skip_ws(s);
    if (**s == '/') {
        (*s)++;
        skip_ws(s);
        long den = strtol(*s, &end, 10);
        if (end == *s || den == 0) return false;
        *s = end;
        rat_set_si(out, num, den);
    } else {
        rat_set_si(out, num, 1);
    }
    return true;
}

Vector* vector_parse(const char* str) {
    const char* s = str;
    skip_ws(&s);
    if (*s != '[') return NULL;
    s++;

    int cap = 8, n = 0;
    Rational* data = malloc(cap * sizeof(Rational));

    while (1) {
        skip_ws(&s);
        if (*s == ']') { s++; break; }
        if (n > 0) {
            if (*s != ',') { for(int i=0;i<n;i++) rat_clear(&data[i]); free(data); return NULL; }
            s++;
        }
        if (n >= cap) { cap *= 2; data = realloc(data, cap * sizeof(Rational)); }
        rat_init(&data[n]);
        if (!parse_rational_lit(&s, &data[n])) {
            for(int i=0;i<=n;i++) rat_clear(&data[i]);
            free(data); return NULL;
        }
        n++;
    }

    Vector* v = vector_create(n);
    for (int i = 0; i < n; i++) {
        rat_set(&v->data[i], &data[i]);
        rat_clear(&data[i]);
    }
    free(data);
    return v;
}

Matrix* matrix_parse(const char* str) {
    const char* s = str;
    skip_ws(&s);
    if (*s != '[') return NULL;
    s++;

    int col_count = -1, rows = 0;
    int total = 0, data_cap = 16;
    Rational* data = malloc(data_cap * sizeof(Rational));

    while (1) {
        skip_ws(&s);
        if (*s == ']') { s++; break; }
        if (rows > 0) { skip_ws(&s); if (*s == ',') s++; }
        skip_ws(&s);
        if (*s != '[') { goto fail; }
        s++;

        int cols_this_row = 0;
        while (1) {
            skip_ws(&s);
            if (*s == ']') { s++; break; }
            if (cols_this_row > 0) { if (*s != ',') goto fail; s++; }
            if (total >= data_cap) { data_cap *= 2; data = realloc(data, data_cap * sizeof(Rational)); }
            rat_init(&data[total]);
            if (!parse_rational_lit(&s, &data[total])) { rat_clear(&data[total]); goto fail; }
            total++;
            cols_this_row++;
        }
        if (col_count == -1) col_count = cols_this_row;
        else if (cols_this_row != col_count) goto fail;
        rows++;
    }

    if (rows == 0 || col_count <= 0) goto fail;

    {
        Matrix* m = matrix_create(rows, col_count);
        for (int i = 0; i < total; i++) {
            rat_set(&m->data[i], &data[i]);
            rat_clear(&data[i]);
        }
        free(data);
        return m;
    }

fail:
    for (int i = 0; i < total; i++) rat_clear(&data[i]);
    free(data);
    return NULL;
}
