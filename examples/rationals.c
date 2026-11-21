// Exact rationals: literals are parsed exactly (0.1 is 1/10), and nothing is rounded.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, const ps_rational* r) {
    char* s = ps_rat_to_string(r);
    printf("%-22s %s\n", label, s);
    free(s);
}

int main(void) {
    ps_rational a, b, sum, big;
    ps_rat_init(&a); ps_rat_init(&b); ps_rat_init(&sum); ps_rat_init(&big);
    ps_rat_parse(&a, "0.1");                  // exactly 1/10
    ps_rat_set_si(&b, 1, 3);
    ps_rat_add(&sum, &a, &b);
    show("0.1 + 1/3 =", &sum);                     // 13/30
    ps_rat_parse(&big, "2.5e30");
    ps_rat_mul(&big, &big, &big);
    show("(2.5e30)^2 =", &big);                    // every digit kept
    char* dec = ps_rat_to_decimal(&b, 20);
    printf("%-22s %s\n", "1/3 to 20 places:", dec);
    free(dec);
    ps_rat_clear(&a); ps_rat_clear(&b); ps_rat_clear(&sum); ps_rat_clear(&big);
    return 0;
}
