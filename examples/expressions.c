// Building expressions in code instead of parsing text, then simplifying, differentiating
// and evaluating them.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, const Expr* e) {
    char* s = expr_to_string(e);
    printf("%-14s %s\n", label, s);
    free(s);
}

int main(void) {
    // f = (x^2 - 1)/(x - 1) + 1/x, built from constructors (each takes its arguments)
    Expr* x2 = create_pow(create_var("x"), create_const(2));
    Expr* f = create_add(create_div(create_sub(x2, create_const(1)), create_sub(create_var("x"), create_const(1))),
                         create_div(create_const(1), create_var("x")));
    show("f =", f);
    Expr* s = simplify(f);                         // a reduced rational function
    show("simplified:", s);
    Expr* d = differentiate(f, "x");
    show("f' =", d);
    bool ok = true;
    printf("%-14s %.12g\n", "f(2.5) =", expr_eval(f, 2.5, &ok));

    // with exp: canonical Laurent form in exp(x) through simplify_exp
    const char* why = NULL;
    Expr* g = parse_input("(exp(2x) - 1)/(exp(x) - 1)");
    Expr* gs = simplify_exp(g, &why);
    show("simplify_exp:", gs);
    free_expr(f); free_expr(s); free_expr(d); free_expr(g); free_expr(gs);
    return 0;
}
