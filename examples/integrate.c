// Indefinite integrals: rational functions (ps_integrate), exp/ln/sin/cos (ps_integrate_exp, the
// Risch algorithm, which also proves when no elementary antiderivative exists), and one square
// root of a quadratic (ps_integrate_sqrt).
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static void rational(const char* text) {
    ps_expr* f = ps_parse(text);
    ps_expr* F = ps_integrate(f, "x");                    // NULL with ps_integrate_error() on failure
    char* s = F ? ps_expr_to_string(F) : NULL;
    printf("∫ %-22s = %s\n", text, s ? s : ps_integrate_error());
    free(s); ps_expr_free(F); ps_expr_free(f);
}

static void risch(const char* text) {
    ps_expr* f = ps_parse(text);
    ps_risch_result r = ps_integrate_exp(f);
    char* s = r.result ? ps_expr_to_string(r.result) : NULL;
    if (r.status == PS_RISCH_OK) printf("∫ %-22s = %s\n", text, s);
    else if (r.status == PS_RISCH_NONELEMENTARY) printf("∫ %-22s : no elementary antiderivative (the part %s)\n", text, s);
    else printf("∫ %-22s : unsupported (%s)\n", text, r.reason);
    free(s); ps_expr_free(r.result); ps_expr_free(f);
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

    const char* roots[] = { "sqrt(x^2+1)", "1/sqrt(1-x^2)", "x/sqrt(x+1)" };
    for (size_t i = 0; i < 3; i++) {              // ps_integrate_sqrt: PS_RISCH_OK or PS_RISCH_UNSUPPORTED
        ps_expr* f = ps_parse(roots[i]);
        ps_risch_result r = ps_integrate_sqrt(f);
        char* s = r.result ? ps_expr_to_string(r.result) : NULL;
        printf("∫ %-22s = %s\n", roots[i], s ? s : r.reason);
        free(s); ps_expr_free(r.result); ps_expr_free(f);
    }

    const char* why = NULL;                       // ps_differentiate_exp: the same expressions
    ps_expr* g = ps_parse("x*exp(x)*sin(x)");
    ps_expr* dg = ps_differentiate_exp(g, &why);
    char* s = ps_expr_to_string(dg);
    printf("d/dx x*exp(x)*sin(x) = %s\n", s);
    free(s); ps_expr_free(dg); ps_expr_free(g);
    return 0;
}
