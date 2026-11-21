// Polynomials over Q: build from text or term by term, then arithmetic, gcd and calculus.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

// text -> polynomial (NULL if the text is not a polynomial)
static ps_poly* poly_of(const char* text) {
    ps_expr* raw = ps_parse(text);
    ps_expr* e = raw ? ps_canonicalize(raw) : NULL;
    ps_poly* p = ps_expr_poly(e) ? ps_poly_copy(ps_expr_poly(e)) : NULL;
    ps_expr_free(raw); ps_expr_free(e);
    return p;
}

static void show(const char* label, const ps_poly* p) {
    char* s = ps_poly_to_string(p);
    printf("%-16s %s\n", label, s);
    free(s);
}

int main(void) {
    ps_poly* a = poly_of("(x - 1)^2 * (x + 2)");
    ps_poly* b = ps_poly_new();          // x^2 - 1, term by term
    ps_rational one = ps_rat_from_int(1), minus_one = ps_rat_from_int(-1);
    ps_poly_add_term(b, &one, 2);
    ps_poly_add_term(b, &minus_one, 0);

    show("a =", a);
    show("b =", b);
    ps_poly* prod = ps_poly_mul(a, b);
    show("a * b =", prod);
    ps_divmod_result qr = ps_poly_divmod(a, b);
    show("a div b =", qr.quotient);
    show("a mod b =", qr.remainder);
    ps_poly* g = ps_poly_gcd(a, b);
    show("gcd(a, b) =", g);
    ps_bezout bz = ps_poly_xgcd(a, b);               // s·a + t·b = gcd
    show("Bezout s =", bz.s);
    show("Bezout t =", bz.t);
    ps_poly* da = ps_poly_derivative(a);
    show("a' =", da);

    ps_rational x0, v;                                // exact evaluation at x = 1/2
    ps_rat_init(&x0); ps_rat_init(&v);
    ps_rat_set_si(&x0, 1, 2);
    ps_poly_eval(&v, a, &x0);
    char* vs = ps_rat_to_string(&v);
    printf("%-16s %s\n", "a(1/2) =", vs);
    free(vs);

    ps_rat_clear(&one); ps_rat_clear(&minus_one); ps_rat_clear(&x0); ps_rat_clear(&v);
    ps_bezout_free(&bz);
    ps_poly_free(a); ps_poly_free(b); ps_poly_free(prod); ps_poly_free(g); ps_poly_free(da);
    ps_poly_free(qr.quotient); ps_poly_free(qr.remainder);
    return 0;
}
