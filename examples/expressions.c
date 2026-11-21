// Building expressions in code instead of parsing text, then simplifying, differentiating
// and evaluating them.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, const ps_expr* e) {
    char* s = ps_expr_to_string(e);
    printf("%-14s %s\n", label, s);
    free(s);
}

int main(void) {
    // f = (x^2 - 1)/(x - 1) + 1/x, built from constructors (each takes its arguments)
    ps_expr* x2 = ps_expr_pow(ps_expr_var("x"), ps_expr_int(2));
    ps_expr* f = ps_expr_add(ps_expr_div(ps_expr_sub(x2, ps_expr_int(1)), ps_expr_sub(ps_expr_var("x"), ps_expr_int(1))),
                         ps_expr_div(ps_expr_int(1), ps_expr_var("x")));
    show("f =", f);
    ps_expr* s = ps_simplify(f);                         // a reduced rational function
    show("simplified:", s);
    ps_expr* d = ps_differentiate(f, "x");
    show("f' =", d);
    bool ok = true;
    printf("%-14s %.12g\n", "f(2.5) =", ps_expr_eval(f, 2.5, &ok));

    // with exp: canonical Laurent form in exp(x) through ps_simplify_exp
    const char* why = NULL;
    ps_expr* g = ps_parse("(exp(2x) - 1)/(exp(x) - 1)");
    ps_expr* gs = ps_simplify_exp(g, &why);
    show("simplify_exp:", gs);
    ps_expr_free(f); ps_expr_free(s); ps_expr_free(d); ps_expr_free(g); ps_expr_free(gs);
    return 0;
}
