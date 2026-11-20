// Indefinite integrals: rational functions (integrate), and exp/ln/sin/cos (integrate_exp,
// the Risch algorithm), which also proves when no elementary antiderivative exists.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static void rational(const char* text) {
    Expr* f = parse_input(text);
    Expr* F = integrate(f, "x");                    // NULL with integrate_error() on failure
    char* s = F ? expr_to_string(F) : NULL;
    printf("∫ %-22s = %s\n", text, s ? s : integrate_error());
    free(s); free_expr(F); free_expr(f);
}

static void risch(const char* text) {
    Expr* f = parse_input(text);
    RischResult r = integrate_exp(f);
    char* s = r.result ? expr_to_string(r.result) : NULL;
    if (r.status == RISCH_OK) printf("∫ %-22s = %s\n", text, s);
    else if (r.status == RISCH_NONELEMENTARY) printf("∫ %-22s : no elementary antiderivative (the part %s)\n", text, s);
    else printf("∫ %-22s : unsupported (%s)\n", text, r.reason);
    free(s); free_expr(r.result); free_expr(f);
}

int main(void) {
    rational("1/((x-1)(x^2+1))");
    rational("x/(x^4+1)");
    rational("1/(x^3-2)");                        // radicals: the roots of x^3 - 2 are known
    rational("1/(x^5+x+3)");                      // no radicals exist: an exact RootSum
    risch("x^2*exp(3x)");
    risch("x^2*ln(x)^3");
    risch("exp(x)/(1+exp(2x))");
    risch("x*exp(x)*sin(x)");
    risch("1/(2+cos(x))");
    risch("exp(x)/(exp(3x)+2)");                  // cubic residues: a RootSum over exp
    risch("exp(x^2)");                            // proved non-elementary
    risch("ln(x)/(x+1)");

    const char* why = NULL;                       // differentiate_exp: the same expressions
    Expr* g = parse_input("x*exp(x)*sin(x)");
    Expr* dg = differentiate_exp(g, &why);
    char* s = expr_to_string(dg);
    printf("d/dx x*exp(x)*sin(x) = %s\n", s);
    free(s); free_expr(dg); free_expr(g);
    return 0;
}
