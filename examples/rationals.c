// Exact rationals: literals are parsed exactly (0.1 is 1/10), and nothing is rounded.
#include "poly_solver.h"
#include <stdio.h>
#include <stdlib.h>

static void show(const char* label, const Rational* r) {
    char* s = rat_to_string(r);
    printf("%-22s %s\n", label, s);
    free(s);
}

int main(void) {
    Rational a, b, sum, big;
    rat_init(&a); rat_init(&b); rat_init(&sum); rat_init(&big);
    rat_parse_literal(&a, "0.1");                  // exactly 1/10
    rat_set_si(&b, 1, 3);
    rat_add(&sum, &a, &b);
    show("0.1 + 1/3 =", &sum);                     // 13/30
    rat_parse_literal(&big, "2.5e30");
    rat_mul(&big, &big, &big);
    show("(2.5e30)^2 =", &big);                    // every digit kept
    char* dec = rat_to_decimal(&b, 20);
    printf("%-22s %s\n", "1/3 to 20 places:", dec);
    free(dec);
    rat_clear(&a); rat_clear(&b); rat_clear(&sum); rat_clear(&big);
    return 0;
}
