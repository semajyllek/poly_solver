// Big numbers: every coefficient is an arbitrary-precision GMP rational, so nothing overflows
// or rounds, in factoring, roots, integrals or matrices.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static ps_expr* simplified(const char* text) {
    ps_expr* raw = ps_parse(text);
    ps_expr* e = ps_canonicalize(raw);
    ps_expr_free(raw);
    return e;
}

int main(void) {
    // (x - 10^30)(x + 2^100): both roots come back exactly
    ps_expr* p = simplified("(x - 1000000000000000000000000000000)(x + 1267650600228229401496703205376)");
    char* ps = ps_expr_to_string(p);
    printf("p        = %s\n", ps);
    ps_factorization* f = ps_factorize(ps_expr_poly(p));
    char* fs = ps_factorization_to_string(f);
    printf("factor p = %s\n", fs);
    ps_solutions* sol = ps_solve(ps_expr_poly(p));
    char* ss = ps_solutions_to_string(sol);
    printf("solve p  = %s\n", ss);

    ps_expr* q = simplified("(x + 1)^60");           // binomial coefficients up to C(60, 30)
    char* qs = ps_expr_to_string(q);
    printf("(x+1)^60 = %.120s...\n", qs);

    ps_expr* r = ps_parse("1/(x^2 + 123456789012345678901234567890)");
    ps_expr* R = ps_integrate(r, "x");
    char* Rs = ps_expr_to_string(R);
    printf("∫ 1/(x^2 + 123456789012345678901234567890) = %s\n", Rs);

    ps_matrix* m = ps_matrix_parse("[[1e40, 3], [7, 1/99999999999999999999]]");
    ps_rational d;
    ps_rat_init(&d);
    ps_matrix_det(&d, m);
    char* ds = ps_rat_to_string(&d);
    printf("det      = %s\n", ds);

    free(ps); free(fs); free(ss); free(qs); free(Rs); free(ds);
    ps_rat_clear(&d);
    ps_matrix_free(m);
    ps_factorization_free(f); ps_solutions_free(sol);
    ps_expr_free(p); ps_expr_free(q); ps_expr_free(r); ps_expr_free(R);
    return 0;
}
