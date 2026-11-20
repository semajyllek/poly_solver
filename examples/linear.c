// Exact linear algebra: matrices of rationals, parsed from bracket notation or filled directly.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, char* s) {
    printf("%-10s %s\n", label, s);
    free(s);
}

int main(void) {
    Matrix* H = matrix_create(3, 3);              // the Hilbert matrix 1/(i + j + 1)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) rat_set_si(&MAT_AT(H, i, j), 1, i + j + 1);
    show("H =", matrix_to_string(H));
    Rational d;
    rat_init(&d);
    matrix_det(&d, H);
    show("det H =", rat_to_string(&d));
    Matrix* Hi = matrix_inverse(H);
    show("H^-1 =", matrix_to_string(Hi));

    Matrix* A = matrix_parse("[[1, 2, 3], [4, 5, 6], [7, 8, 10]]");
    Vector* b = vector_parse("[1, 0.5, 1/3]");
    LinSysResult* x = linsys_solve(A, b);
    show("A x = b:", linsys_result_to_string(x));
    Matrix* S = matrix_parse("[[1, 2, 3], [2, 4, 6], [1, 0, 1]]");
    printf("%-10s %d\n", "rank S =", matrix_rank(S));
    Matrix* R = matrix_rref(S);
    show("rref S =", matrix_to_string(R));

    rat_clear(&d);
    matrix_free(H); matrix_free(Hi); matrix_free(A); matrix_free(S); matrix_free(R);
    vector_free(b); linsys_result_free(x);
    return 0;
}
