#ifndef MATRIX_H
#define MATRIX_H

#include "rational_num.h"
#include <stdbool.h>

typedef ps_matrix Matrix;

typedef ps_vector Vector;

typedef ps_linsys_status LinSysStatus;
#define LINSYS_UNIQUE PS_LINSYS_UNIQUE
#define LINSYS_INCONSISTENT PS_LINSYS_INCONSISTENT
#define LINSYS_INFINITE PS_LINSYS_INFINITE

typedef ps_linsys_result LinSysResult;

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
