// Big numbers: every coefficient is an arbitrary-precision GMP rational, so nothing overflows
// or rounds, in factoring, roots, integrals or matrices.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static Expr* simplified(const char* text) {
    Expr* raw = parse_input(text);
    Expr* e = canonicalize(raw);
    free_expr(raw);
    return e;
}

int main(void) {
    // (x - 10^30)(x + 2^100): both roots come back exactly
    Expr* p = simplified("(x - 1000000000000000000000000000000)(x + 1267650600228229401496703205376)");
    char* ps = expr_to_string(p);
    printf("p        = %s\n", ps);
    Factorization* f = factorize(p->poly);
    char* fs = factorization_to_string(f);
    printf("factor p = %s\n", fs);
    Solutions* sol = solve_polynomial(p->poly);
    char* ss = solutions_to_string(sol);
    printf("solve p  = %s\n", ss);

    Expr* q = simplified("(x + 1)^60");           // binomial coefficients up to C(60, 30)
    char* qs = expr_to_string(q);
    printf("(x+1)^60 = %.120s...\n", qs);

    Expr* r = parse_input("1/(x^2 + 123456789012345678901234567890)");
    Expr* R = integrate(r, "x");
    char* Rs = expr_to_string(R);
    printf("∫ 1/(x^2 + 123456789012345678901234567890) = %s\n", Rs);

    Matrix* m = matrix_parse("[[1e40, 3], [7, 1/99999999999999999999]]");
    Rational d;
    rat_init(&d);
    matrix_det(&d, m);
    char* ds = rat_to_string(&d);
    printf("det      = %s\n", ds);

    free(ps); free(fs); free(ss); free(qs); free(Rs); free(ds);
    rat_clear(&d);
    matrix_free(m);
    free_factorization(f); free_solutions(sol);
    free_expr(p); free_expr(q); free_expr(r); free_expr(R);
    return 0;
}
