#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include "expression.h"
#include "polynomial.h"
#include "factor.h"
#include "solver.h"

static void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s <action> \"<expression>\"\n", program_name);
    fprintf(stderr, "Actions:\n");
    fprintf(stderr, "  simplify      - Simplify the expression\n");
    fprintf(stderr, "  factor        - Factor the expression\n");
    fprintf(stderr, "  differentiate - Differentiate the expression\n");
    fprintf(stderr, "  integrate     - Integrate the expression\n");
    fprintf(stderr, "  solve         - Solve the equation (= 0)\n");
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        print_usage(argv[0]);
        return 1;
    }

    const char* action = argv[1];
    const char* input = argv[2];

    // parse and canonicalize the input
    Expr* expr = parse_and_simplify(input);
    if (expr == NULL) {
        fprintf(stderr, "Error: could not parse expression.\n");
        return 1;
    }

    if (strcmp(action, "simplify") == 0) {
        char* s = expr_to_string(expr);
        printf("%s\n", s);
        free(s);

    } else if (strcmp(action, "factor") == 0) {
        if (expr->type != EXPR_POLY) {
            fprintf(stderr, "Error: can only factor polynomial expressions.\n");
            free_expr(expr);
            return 1;
        }
        Factorization* f = factorize(expr->poly);
        char* s = factorization_to_string(f);
        printf("%s\n", s);
        free(s);
        free_factorization(f);

    } else if (strcmp(action, "differentiate") == 0) {
        Expr* deriv = differentiate(expr, "x");
        if (!deriv) {
            fprintf(stderr, "Error: could not differentiate expression.\n");
            free_expr(expr);
            return 1;
        }
        char* s = expr_to_string(deriv);
        printf("%s\n", s);
        free(s);
        free_expr(deriv);

    } else if (strcmp(action, "integrate") == 0) {
        Expr* integral = integrate(expr, "x");
        if (!integral) {
            fprintf(stderr, "Error: could not integrate expression.\n");
            free_expr(expr);
            return 1;
        }
        char* s = expr_to_string(integral);
        printf("%s\n", s);
        free(s);
        free_expr(integral);

    } else if (strcmp(action, "solve") == 0) {
        if (expr->type != EXPR_POLY) {
            fprintf(stderr, "Error: can only solve polynomial equations.\n");
            free_expr(expr);
            return 1;
        }
        Solutions* sol = solve_polynomial(expr->poly);
        char* s = solutions_to_string(sol);
        printf("%s\n", s);
        free(s);
        free_solutions(sol);

    } else {
        fprintf(stderr, "Error: unknown action '%s'.\n", action);
        free_expr(expr);
        return 1;
    }

    free_expr(expr);
    return 0;
}
