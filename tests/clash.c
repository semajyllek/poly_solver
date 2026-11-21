// A program that defines names the library uses internally (Expr, integrate, poly_add, RISCH_OK,
// create_add) and still links and runs against the library: only ps_ names are shared.
#include "polysolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Expr { int mine; } Expr;          // the library's internal type name
#define RISCH_OK 42                               // and an internal enum constant
static int calls;
int integrate(int x) { return calls += x; }       // internal function names, defined here
int poly_add(int a, int b) { return a + b; }
void* create_add(void* a, void* b) { (void)b; return a; }
int free_expr = 7;                                // a global variable with an internal name

int main(void) {
    Expr e = { RISCH_OK };
    ps_expr* f = ps_parse_simplify("1/(x^2+1)");
    ps_expr* F = ps_integrate(f, "x");                 // integrate leaves f alone
    char* s = ps_expr_to_string(F);
    int ok = strcmp(s, "atan(x)") == 0 && integrate(1) == 1 && poly_add(2, 3) == 5 && e.mine == 42
             && create_add(&e, NULL) == &e && free_expr == 7 && strcmp(ps_version(), PS_VERSION) == 0;
    printf("clash: %s (%s, version %s)\n", ok ? "ok" : "FAILED", s, ps_version());
    free(s);
    ps_expr_free(f);
    ps_expr_free(F);
    return ok ? 0 : 1;
}
