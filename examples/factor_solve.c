// Factoring over Q, and solving: exact roots where they exist, certified intervals elsewhere.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static Polynomial* poly_of(const char* text) {
    Expr* raw = parse_input(text);
    Expr* e = raw ? canonicalize(raw) : NULL;
    Polynomial* p = e && e->type == EXPR_POLY ? poly_copy(e->poly) : NULL;
    free_expr(raw); free_expr(e);
    return p;
}

int main(void) {
    Polynomial* p = poly_of("2x^6 - 2x^4 + 8x^2 - 8");
    Factorization* f = factorize(p);               // p == content · Π factors^mult
    char* fs = factorization_to_string(f);
    printf("factor: %s\n", fs);
    free(fs);
    for (int i = 0; i < f->count; i++) {
        char* s = poly_to_string(f->factors[i]);
        printf("  factor %d: %s (multiplicity %d)\n", i + 1, s, f->multiplicities[i]);
        free(s);
    }

    Solutions* sol = solve_polynomial(p);
    char* ss = solutions_to_string(sol);
    printf("solve:  %s\n", ss);
    free(ss);

    Polynomial* q = poly_of("x^5 - 3x + 1");     // no radicals: certified real-root intervals
    RootIntervalList L = isolate_real_roots(q);
    Rational w;
    rat_init(&w);
    rat_set_si(&w, 1, 1000000);
    for (int i = 0; i < L.count; i++) {
        refine_root(&L.items[i], q, &w);          // shrink to width 1e-6, exactly
        char *lo = rat_to_decimal(&L.items[i].lo, 8), *hi = rat_to_decimal(&L.items[i].hi, 8);
        printf("  root of x^5 - 3x + 1 in (%s, %s)\n", lo, hi);
        free(lo); free(hi);
    }
    rat_clear(&w);
    root_interval_list_free(&L);
    free_solutions(sol); free_factorization(f);
    free_polynomial(p); free_polynomial(q);
    return 0;
}
