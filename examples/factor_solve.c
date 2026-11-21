// Factoring over Q, and solving: exact roots where they exist, certified intervals elsewhere.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static ps_poly* poly_of(const char* text) {
    ps_expr* raw = ps_parse(text);
    ps_expr* e = raw ? ps_canonicalize(raw) : NULL;
    ps_poly* p = ps_expr_poly(e) ? ps_poly_copy(ps_expr_poly(e)) : NULL;
    ps_expr_free(raw); ps_expr_free(e);
    return p;
}

int main(void) {
    ps_poly* p = poly_of("2x^6 - 2x^4 + 8x^2 - 8");
    ps_factorization* f = ps_factorize(p);               // p == content · Π factors^mult
    char* fs = ps_factorization_to_string(f);
    printf("factor: %s\n", fs);
    free(fs);
    for (int i = 0; i < f->count; i++) {
        char* s = ps_poly_to_string(f->factors[i]);
        printf("  factor %d: %s (multiplicity %d)\n", i + 1, s, f->multiplicities[i]);
        free(s);
    }

    ps_solutions* sol = ps_solve(p);
    char* ss = ps_solutions_to_string(sol);
    printf("solve:  %s\n", ss);
    free(ss);

    ps_poly* q = poly_of("x^5 - 3x + 1");     // no radicals: certified real-root intervals
    ps_root_list L = ps_real_roots(q);
    ps_rational w;
    ps_rat_init(&w);
    ps_rat_set_si(&w, 1, 1000000);
    for (int i = 0; i < L.count; i++) {
        ps_root_refine(&L.items[i], q, &w);          // shrink to width 1e-6, exactly
        char *lo = ps_rat_to_decimal(&L.items[i].lo, 8), *hi = ps_rat_to_decimal(&L.items[i].hi, 8);
        printf("  root of x^5 - 3x + 1 in (%s, %s)\n", lo, hi);
        free(lo); free(hi);
    }
    ps_rat_clear(&w);
    ps_root_list_free(&L);
    ps_solutions_free(sol); ps_factorization_free(f);
    ps_poly_free(p); ps_poly_free(q);
    return 0;
}
