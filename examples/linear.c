// Exact linear algebra: matrices of rationals, parsed from bracket notation or filled directly.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, char* s) {
    printf("%-10s %s\n", label, s);
    free(s);
}

int main(void) {
    ps_matrix* H = ps_matrix_new(3, 3);              // the Hilbert matrix 1/(i + j + 1)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) ps_rat_set_si(&PS_MAT_AT(H, i, j), 1, i + j + 1);
    show("H =", ps_matrix_to_string(H));
    ps_rational d;
    ps_rat_init(&d);
    ps_matrix_det(&d, H);
    show("det H =", ps_rat_to_string(&d));
    ps_matrix* Hi = ps_matrix_inverse(H);
    show("H^-1 =", ps_matrix_to_string(Hi));

    ps_matrix* A = ps_matrix_parse("[[1, 2, 3], [4, 5, 6], [7, 8, 10]]");
    ps_vector* b = ps_vector_parse("[1, 0.5, 1/3]");
    ps_linsys_result* x = ps_linsys_solve(A, b);
    show("A x = b:", ps_linsys_to_string(x));
    ps_matrix* S = ps_matrix_parse("[[1, 2, 3], [2, 4, 6], [1, 0, 1]]");
    printf("%-10s %d\n", "rank S =", ps_matrix_rank(S));
    ps_matrix* R = ps_matrix_rref(S);
    show("rref S =", ps_matrix_to_string(R));

    ps_rat_clear(&d);
    ps_matrix_free(H); ps_matrix_free(Hi); ps_matrix_free(A); ps_matrix_free(S); ps_matrix_free(R);
    ps_vector_free(b); ps_linsys_free(x);
    return 0;
}
