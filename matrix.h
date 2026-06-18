#ifndef MATRIX_H
#define MATRIX_H

#include "rational_num.h"
#include <stdbool.h>

typedef struct {
    Rational* data; // row-major: element (r,c) at data[r * cols + c]
    int rows;
    int cols;
} Matrix;

typedef struct {
    Rational* data;
    int n;
} Vector;

typedef enum {
    LINSYS_UNIQUE,
    LINSYS_INCONSISTENT,
    LINSYS_INFINITE
} LinSysStatus;

typedef struct {
    LinSysStatus status;
    Vector* solution; // NULL if inconsistent
    int rank;
} LinSysResult;

#define MAT_AT(m, r, c) ((m)->data[(r) * (m)->cols + (c)])

// construction / destruction
Matrix* matrix_create(int rows, int cols);
Matrix* matrix_copy(const Matrix* m);
void    matrix_free(Matrix* m);

Vector* vector_create(int n);
Vector* vector_copy(const Vector* v);
void    vector_free(Vector* v);

// core operations
Matrix* matrix_rref(const Matrix* m);
int     matrix_rank(const Matrix* m);
void matrix_det(Rational* result, const Matrix* m);
Matrix* matrix_inverse(const Matrix* m);
LinSysResult* linsys_solve(const Matrix* A, const Vector* b);

// output
char* matrix_to_string(const Matrix* m);
char* vector_to_string(const Vector* v);
char* linsys_result_to_string(const LinSysResult* r);
void  linsys_result_free(LinSysResult* r);

// CLI parsing
Matrix* matrix_parse(const char* str);
Vector* vector_parse(const char* str);

#endif
