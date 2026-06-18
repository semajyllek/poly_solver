#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

// --- construction / destruction ---

Matrix* matrix_create(int rows, int cols) {
    Matrix* m = malloc(sizeof(Matrix));
    m->rows = rows;
    m->cols = cols;
    m->data = calloc(rows * cols, sizeof(Rational));
    // calloc zeros memory; Rational {0, 0} isn't normalized, so init properly
    for (int i = 0; i < rows * cols; i++)
        m->data[i] = rat_zero();
    return m;
}

Matrix* matrix_copy(const Matrix* m) {
    Matrix* copy = matrix_create(m->rows, m->cols);
    memcpy(copy->data, m->data, m->rows * m->cols * sizeof(Rational));
    return copy;
}

void matrix_free(Matrix* m) {
    if (!m) return;
    free(m->data);
    free(m);
}

Vector* vector_create(int n) {
    Vector* v = malloc(sizeof(Vector));
    v->n = n;
    v->data = calloc(n, sizeof(Rational));
    for (int i = 0; i < n; i++)
        v->data[i] = rat_zero();
    return v;
}

Vector* vector_copy(const Vector* v) {
    Vector* copy = vector_create(v->n);
    memcpy(copy->data, v->data, v->n * sizeof(Rational));
    return copy;
}

void vector_free(Vector* v) {
    if (!v) return;
    free(v->data);
    free(v);
}

// --- internal: Gaussian elimination to RREF ---
// operates on a mutable matrix, returns number of row swaps
static int gauss_eliminate(Matrix* m, int augmented_cols) {
    int pivot_row = 0;
    int pivot_col = 0;
    int swaps = 0;
    int coeff_cols = m->cols - augmented_cols;

    while (pivot_row < m->rows && pivot_col < coeff_cols) {
        // find nonzero entry in this column at or below pivot_row
        int found = -1;
        for (int r = pivot_row; r < m->rows; r++) {
            if (!rat_is_zero(MAT_AT(m, r, pivot_col))) {
                found = r;
                break;
            }
        }

        if (found == -1) {
            pivot_col++;
            continue;
        }

        // swap rows if needed
        if (found != pivot_row) {
            for (int c = 0; c < m->cols; c++) {
                Rational tmp = MAT_AT(m, pivot_row, c);
                MAT_AT(m, pivot_row, c) = MAT_AT(m, found, c);
                MAT_AT(m, found, c) = tmp;
            }
            swaps++;
        }

        // scale pivot row so pivot = 1
        Rational pivot_val = MAT_AT(m, pivot_row, pivot_col);
        for (int c = pivot_col; c < m->cols; c++)
            MAT_AT(m, pivot_row, c) = rat_div(MAT_AT(m, pivot_row, c), pivot_val);

        // eliminate all other rows
        for (int r = 0; r < m->rows; r++) {
            if (r == pivot_row) continue;
            Rational factor = MAT_AT(m, r, pivot_col);
            if (rat_is_zero(factor)) continue;
            for (int c = pivot_col; c < m->cols; c++)
                MAT_AT(m, r, c) = rat_sub(MAT_AT(m, r, c), rat_mul(factor, MAT_AT(m, pivot_row, c)));
        }

        pivot_row++;
        pivot_col++;
    }

    return swaps;
}

// --- core operations ---

Matrix* matrix_rref(const Matrix* m) {
    Matrix* result = matrix_copy(m);
    gauss_eliminate(result, 0);
    return result;
}

int matrix_rank(const Matrix* m) {
    Matrix* rref = matrix_rref(m);
    int rank = 0;
    for (int r = 0; r < rref->rows; r++) {
        bool nonzero = false;
        for (int c = 0; c < rref->cols; c++) {
            if (!rat_is_zero(MAT_AT(rref, r, c))) {
                nonzero = true;
                break;
            }
        }
        if (nonzero) rank++;
    }
    matrix_free(rref);
    return rank;
}

Rational matrix_det(const Matrix* m) {
    assert(m->rows == m->cols);
    int n = m->rows;
    Matrix* work = matrix_copy(m);
    Rational det = rat_one();
    int swaps = 0;

    for (int col = 0; col < n; col++) {
        // find nonzero pivot
        int found = -1;
        for (int r = col; r < n; r++) {
            if (!rat_is_zero(MAT_AT(work, r, col))) {
                found = r;
                break;
            }
        }
        if (found == -1) {
            matrix_free(work);
            return rat_zero();
        }

        if (found != col) {
            for (int c = 0; c < n; c++) {
                Rational tmp = MAT_AT(work, col, c);
                MAT_AT(work, col, c) = MAT_AT(work, found, c);
                MAT_AT(work, found, c) = tmp;
            }
            swaps++;
        }

        Rational pivot = MAT_AT(work, col, col);
        det = rat_mul(det, pivot);

        // forward eliminate below pivot only
        for (int r = col + 1; r < n; r++) {
            Rational factor = rat_div(MAT_AT(work, r, col), pivot);
            for (int c = col; c < n; c++)
                MAT_AT(work, r, c) = rat_sub(MAT_AT(work, r, c), rat_mul(factor, MAT_AT(work, col, c)));
        }
    }

    matrix_free(work);
    if (swaps % 2 == 1) det = rat_neg(det);
    return det;
}

Matrix* matrix_inverse(const Matrix* m) {
    assert(m->rows == m->cols);
    int n = m->rows;

    // build augmented [A | I]
    Matrix* aug = matrix_create(n, 2 * n);
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++)
            MAT_AT(aug, r, c) = MAT_AT(m, r, c);
        MAT_AT(aug, r, n + r) = rat_one();
    }

    gauss_eliminate(aug, n);

    // check left half is identity
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            Rational expected = (r == c) ? rat_one() : rat_zero();
            if (!rat_eq(MAT_AT(aug, r, c), expected)) {
                matrix_free(aug);
                return NULL;
            }
        }
    }

    // extract right half
    Matrix* inv = matrix_create(n, n);
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
            MAT_AT(inv, r, c) = MAT_AT(aug, r, n + c);

    matrix_free(aug);
    return inv;
}

LinSysResult* linsys_solve(const Matrix* A, const Vector* b) {
    assert(A->rows == b->n);

    LinSysResult* res = malloc(sizeof(LinSysResult));
    res->solution = NULL;
    res->rank = 0;

    // build augmented [A | b]
    Matrix* aug = matrix_create(A->rows, A->cols + 1);
    for (int r = 0; r < A->rows; r++) {
        for (int c = 0; c < A->cols; c++)
            MAT_AT(aug, r, c) = MAT_AT(A, r, c);
        MAT_AT(aug, r, A->cols) = b->data[r];
    }

    gauss_eliminate(aug, 1);

    // count rank (pivot rows in coefficient columns)
    int rank = 0;
    for (int r = 0; r < aug->rows; r++) {
        bool has_nonzero = false;
        for (int c = 0; c < A->cols; c++) {
            if (!rat_is_zero(MAT_AT(aug, r, c))) {
                has_nonzero = true;
                break;
            }
        }
        if (has_nonzero) rank++;
    }
    res->rank = rank;

    // check for inconsistency: row [0 0 ... 0 | c] with c != 0
    for (int r = 0; r < aug->rows; r++) {
        bool all_zero_coeffs = true;
        for (int c = 0; c < A->cols; c++) {
            if (!rat_is_zero(MAT_AT(aug, r, c))) {
                all_zero_coeffs = false;
                break;
            }
        }
        if (all_zero_coeffs && !rat_is_zero(MAT_AT(aug, r, A->cols))) {
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

    // unique solution: read from RREF
    res->status = LINSYS_UNIQUE;
    res->solution = vector_create(A->cols);
    for (int i = 0; i < A->cols; i++)
        res->solution->data[i] = MAT_AT(aug, i, A->cols);

    matrix_free(aug);
    return res;
}

// --- output ---

char* matrix_to_string(const Matrix* m) {
    size_t cap = 128, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';

    for (int r = 0; r < m->rows; r++) {
        const char* row_start = (r == 0) ? "[[" : " [";
        while (len + strlen(row_start) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, row_start);
        len += strlen(row_start);

        for (int c = 0; c < m->cols; c++) {
            char* rs = rat_to_string(MAT_AT(m, r, c));
            const char* sep = (c < m->cols - 1) ? ", " : "";
            size_t needed = strlen(rs) + strlen(sep) + 1;
            while (len + needed > cap) { cap *= 2; buf = realloc(buf, cap); }
            strcat(buf, rs);
            strcat(buf, sep);
            len += strlen(rs) + strlen(sep);
            free(rs);
        }

        const char* row_end = (r < m->rows - 1) ? "]\n" : "]]";
        while (len + strlen(row_end) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, row_end);
        len += strlen(row_end);
    }

    return buf;
}

char* vector_to_string(const Vector* v) {
    size_t cap = 64, len = 0;
    char* buf = malloc(cap);
    buf[0] = '\0';
    strcat(buf, "["); len++;

    for (int i = 0; i < v->n; i++) {
        char* rs = rat_to_string(v->data[i]);
        const char* sep = (i < v->n - 1) ? ", " : "";
        size_t needed = strlen(rs) + strlen(sep) + 1;
        while (len + needed > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, rs);
        strcat(buf, sep);
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
        snprintf(buf, cap, "Infinitely many solutions (rank = %d, %d variables)",
                 r->rank, r->rank); // simplified
        return buf;
    }

    // unique solution
    size_t len = 0;
    buf[0] = '\0';
    for (int i = 0; i < r->solution->n; i++) {
        char* rs = rat_to_string(r->solution->data[i]);
        char entry[128];
        snprintf(entry, sizeof(entry), "x%d = %s", i + 1, rs);
        free(rs);
        size_t elen = strlen(entry);
        const char* sep = (i < r->solution->n - 1) ? ", " : "";
        while (len + elen + strlen(sep) + 1 > cap) { cap *= 2; buf = realloc(buf, cap); }
        strcat(buf, entry);
        strcat(buf, sep);
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

// skip whitespace
static void skip_ws(const char** s) {
    while (isspace(**s)) (*s)++;
}

// parse a rational number: integer or p/q
static bool parse_rational(const char** s, Rational* out) {
    skip_ws(s);
    char* end;
    long long num = strtoll(*s, &end, 10);
    if (end == *s) return false;
    *s = end;
    skip_ws(s);
    if (**s == '/') {
        (*s)++;
        skip_ws(s);
        long long den = strtoll(*s, &end, 10);
        if (end == *s || den == 0) return false;
        *s = end;
        *out = rat_create(num, den);
    } else {
        *out = rat_from_int(num);
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
            if (*s != ',') { free(data); return NULL; }
            s++;
        }
        Rational val;
        if (!parse_rational(&s, &val)) { free(data); return NULL; }
        if (n >= cap) { cap *= 2; data = realloc(data, cap * sizeof(Rational)); }
        data[n++] = val;
    }

    Vector* v = vector_create(n);
    memcpy(v->data, data, n * sizeof(Rational));
    free(data);
    return v;
}

Matrix* matrix_parse(const char* str) {
    const char* s = str;
    skip_ws(&s);
    if (*s != '[') return NULL;
    s++;

    int col_count = -1, rows = 0;
    Rational* data = NULL;
    int total = 0, data_cap = 16;
    data = malloc(data_cap * sizeof(Rational));

    while (1) {
        skip_ws(&s);
        if (*s == ']') { s++; break; }
        if (rows > 0) {
            skip_ws(&s);
            if (*s == ',') s++;
        }

        // parse one row: [a, b, c]
        skip_ws(&s);
        if (*s != '[') { free(data); return NULL; }
        s++;

        int cols_this_row = 0;
        while (1) {
            skip_ws(&s);
            if (*s == ']') { s++; break; }
            if (cols_this_row > 0) {
                if (*s != ',') { free(data); return NULL; }
                s++;
            }
            Rational val;
            if (!parse_rational(&s, &val)) { free(data); return NULL; }
            if (total >= data_cap) { data_cap *= 2; data = realloc(data, data_cap * sizeof(Rational)); }
            data[total++] = val;
            cols_this_row++;
        }

        if (col_count == -1) col_count = cols_this_row;
        else if (cols_this_row != col_count) { free(data); return NULL; } // ragged
        rows++;
    }

    if (rows == 0 || col_count <= 0) { free(data); return NULL; }

    Matrix* m = matrix_create(rows, col_count);
    memcpy(m->data, data, rows * col_count * sizeof(Rational));
    free(data);
    return m;
}
