// Definite integrals: exact values (constants folded) checked against quadrature at run time.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static void definite(const char* text, const char* a_text, const char* b_text) {
    Rational a, b;
    rat_init(&a); rat_init(&b);
    Bounds ab = { { &a, 0, false }, { &b, 0, false } };
    Expr* f = parse_integrand(text);
    parse_endpoint(&a, &ab.a, a_text);            // "1/2", "pi/4", "inf", "-inf"
    parse_endpoint(&b, &ab.b, b_text);
    DefiniteResult r = integrate_definite(f, &ab);
    printf("∫ %-15s on [%s, %s]: ", text, a_text, b_text);
    if (r.error) printf("error: %s\n", r.error);
    else if (r.exact) { char* s = expr_to_string(r.exact); printf("= %s ~ %.12g\n", s, r.approx); free(s); }
    else printf("~ %.12g (%s)\n", r.approx, r.note);
    definite_free(&r);
    free_expr(f);
    rat_clear(&a); rat_clear(&b);
}

int main(void) {
    definite("1/(x^2+x+1)", "0", "1");
    definite("1/(x^4+1)", "-inf", "inf");
    definite("x^2*sin(x)", "0", "pi");
    definite("exp(-x)*cos(x)", "0", "inf");
    definite("1/(2+cos(x))", "0", "2pi");          // the jump of atan(tan(x/2)/√3) is corrected
    definite("exp(-x^2)", "-inf", "inf");          // no elementary antiderivative: the number
    definite("1/x", "-1", "1");                    // a pole in range
    return 0;
}
