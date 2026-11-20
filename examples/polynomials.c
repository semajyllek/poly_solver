// Polynomials over Q: build from text or term by term, then arithmetic, gcd and calculus.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

// text -> polynomial (NULL if the text is not a polynomial)
static Polynomial* poly_of(const char* text) {
    Expr* raw = parse_input(text);
    Expr* e = raw ? canonicalize(raw) : NULL;
    Polynomial* p = e && e->type == EXPR_POLY ? poly_copy(e->poly) : NULL;
    free_expr(raw); free_expr(e);
    return p;
}

static void show(const char* label, const Polynomial* p) {
    char* s = poly_to_string(p);
    printf("%-16s %s\n", label, s);
    free(s);
}

int main(void) {
    Polynomial* a = poly_of("(x - 1)^2 * (x + 2)");
    Polynomial* b = create_polynomial();          // x^2 - 1, term by term
    Rational one = rat_from_int(1), minus_one = rat_from_int(-1);
    add_term(b, &one, 2);
    add_term(b, &minus_one, 0);

    show("a =", a);
    show("b =", b);
    Polynomial* prod = poly_multiply(a, b);
    show("a * b =", prod);
    PolyDivResult qr = poly_divmod(a, b);
    show("a div b =", qr.quotient);
    show("a mod b =", qr.remainder);
    Polynomial* g = poly_gcd(a, b);
    show("gcd(a, b) =", g);
    PolyBezout bz = poly_xgcd(a, b);               // s·a + t·b = gcd
    show("Bezout s =", bz.s);
    show("Bezout t =", bz.t);
    Polynomial* da = poly_derivative(a);
    show("a' =", da);

    Rational x0, v;                                // exact evaluation at x = 1/2
    rat_init(&x0); rat_init(&v);
    rat_set_si(&x0, 1, 2);
    poly_eval_rat(&v, a, &x0);
    char* vs = rat_to_string(&v);
    printf("%-16s %s\n", "a(1/2) =", vs);
    free(vs);

    rat_clear(&one); rat_clear(&minus_one); rat_clear(&x0); rat_clear(&v);
    free_bezout(&bz);
    free_polynomial(a); free_polynomial(b); free_polynomial(prod); free_polynomial(g); free_polynomial(da);
    free_polynomial(qr.quotient); free_polynomial(qr.remainder);
    return 0;
}
