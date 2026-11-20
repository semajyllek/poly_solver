#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include "expression.h"
#include "polynomial.h"
#include "factor.h"
#include "solver.h"
#include "matrix.h"
#include "definite.h"
#include "risch.h"

static void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s <action> \"<expression>\"\n", program_name);
    fprintf(stderr, "Actions:\n");
    fprintf(stderr, "  simplify      - Simplify the expression\n");
    fprintf(stderr, "  factor        - Factor the expression\n");
    fprintf(stderr, "  differentiate - Differentiate the expression\n");
    fprintf(stderr, "  integrate     - Integrate the expression (add bounds a b for a definite integral; inf/-inf allowed)\n");
    fprintf(stderr, "  solve         - Solve polynomial equation (= 0)\n");
    fprintf(stderr, "  solve_system  - Solve linear system: solve_system \"[[A]]\" \"[b]\"\n");
    fprintf(stderr, "  det           - Determinant of a matrix\n");
    fprintf(stderr, "  inverse       - Inverse of a matrix\n");
    fprintf(stderr, "  rank          - Rank of a matrix\n");
    fprintf(stderr, "  rref          - Reduced row echelon form\n");
}

// matrix actions that take a single matrix argument
static int handle_matrix_action(const char* action, const char* input) {
    Matrix* m = matrix_parse(input);
    if (!m) {
        fprintf(stderr, "Error: could not parse matrix.\n");
        return 1;
    }

    if (strcmp(action, "det") == 0) {
        if (m->rows != m->cols) {
            fprintf(stderr, "Error: determinant requires a square matrix.\n");
            matrix_free(m);
            return 1;
        }
        Rational d;
        rat_init(&d);
        matrix_det(&d, m);
        char* s = rat_to_string(&d);
        printf("%s\n", s);
        free(s);
        rat_clear(&d);
    } else if (strcmp(action, "inverse") == 0) {
        if (m->rows != m->cols) {
            fprintf(stderr, "Error: inverse requires a square matrix.\n");
            matrix_free(m);
            return 1;
        }
        Matrix* inv = matrix_inverse(m);
        if (!inv) {
            printf("Matrix is singular (no inverse).\n");
        } else {
            char* s = matrix_to_string(inv);
            printf("%s\n", s);
            free(s);
            matrix_free(inv);
        }
    } else if (strcmp(action, "rank") == 0) {
        printf("%d\n", matrix_rank(m));
    } else if (strcmp(action, "rref") == 0) {
        Matrix* r = matrix_rref(m);
        char* s = matrix_to_string(r);
        printf("%s\n", s);
        free(s);
        matrix_free(r);
    }

    matrix_free(m);
    return 0;
}

// "exact ~ decimal", or "~ decimal (reason)" when only the number is available
static int handle_definite(const char* f_text, const char* a_text, const char* b_text) {
    Rational a, b;
    rat_init(&a); rat_init(&b);
    Expr* f = parse_integrand(f_text);
    int rc = 1;
    Bounds ab = { { &a, 0, false }, { &b, 0, false } };
    if (!f) fprintf(stderr, "Error: not a polynomial, rational function, or expression in exp, sin and cos.\n");
    else if (!parse_endpoint(&a, &ab.a, a_text) || !parse_endpoint(&b, &ab.b, b_text))
        fprintf(stderr, "Error: bounds must be rational numbers, rational multiples of pi, inf or -inf.\n");
    else {
        DefiniteResult r = integrate_definite(f, &ab);
        if (r.error) fprintf(stderr, "Error: %s.\n", r.error);
        else {
            if (r.exact && r.note) { char* s = expr_to_string(r.exact); printf("%s (%s)\n", s, r.note); free(s); }
            else if (r.exact) { char* s = expr_to_string(r.exact); printf("%s ~ %.10g\n", s, r.approx); free(s); }
            else printf("~ %.10g (%s)\n", r.approx, r.note);
            rc = 0;
        }
        definite_free(&r);
    }
    if (f) free_expr(f);
    rat_clear(&a); rat_clear(&b);
    return rc;
}

// print e and free it
static int print_and_free(Expr* e) {
    char* s = expr_to_string(e);
    printf("%s\n", s);
    free(s);
    free_expr(e);
    return 0;
}

// actions on an expression with exp(...): simplify, differentiate, integrate (Risch)
static int handle_exp(const char* action, Expr* raw) {
    const char* why = NULL;
    Expr* out = NULL;
    int rc = 1;
    if (strcmp(action, "integrate") == 0) {
        RischResult r = integrate_exp(raw);
        if (r.status == RISCH_OK) rc = print_and_free(r.result);
        else if (r.status == RISCH_NONELEMENTARY) {        // a proof, not an error
            char* t = expr_to_string(r.result);
            printf("no elementary antiderivative: the part %s %s\n", t, r.reason + strlen("has no elementary antiderivative "));
            free(t);
            free_expr(r.result);
            rc = 0;
        }
        else fprintf(stderr, "Error: could not integrate: %s.\n", r.reason);
    } else if (strcmp(action, "differentiate") == 0 || strcmp(action, "simplify") == 0) {
        out = action[0] == 'd' ? differentiate_exp(raw, &why) : simplify_exp(raw, &why);
        if (out) rc = print_and_free(out);
        else fprintf(stderr, "Error: %s.\n", why);
    } else {
        fprintf(stderr, "Error: '%s' works on polynomials; exp, ln, sin and cos are supported by simplify, differentiate and integrate.\n", action);
    }
    free_expr(raw);
    return rc;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char* action = argv[1];

    // solve_system needs 4 args: program action matrix vector
    if (strcmp(action, "solve_system") == 0) {
        if (argc != 4) {
            fprintf(stderr, "Usage: %s solve_system \"[[a,b],[c,d]]\" \"[e,f]\"\n", argv[0]);
            return 1;
        }
        Matrix* A = matrix_parse(argv[2]);
        Vector* b = vector_parse(argv[3]);
        if (!A || !b) {
            fprintf(stderr, "Error: could not parse matrix/vector.\n");
            matrix_free(A);
            vector_free(b);
            return 1;
        }
        if (A->rows != b->n) {
            fprintf(stderr, "Error: matrix rows (%d) != vector length (%d).\n", A->rows, b->n);
            matrix_free(A);
            vector_free(b);
            return 1;
        }
        LinSysResult* res = linsys_solve(A, b);
        char* s = linsys_result_to_string(res);
        printf("%s\n", s);
        free(s);
        linsys_result_free(res);
        matrix_free(A);
        vector_free(b);
        return 0;
    }

    // matrix actions that take a single matrix arg
    if (strcmp(action, "det") == 0 || strcmp(action, "inverse") == 0 ||
        strcmp(action, "rank") == 0 || strcmp(action, "rref") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Usage: %s %s \"[[a,b],[c,d]]\"\n", argv[0], action);
            return 1;
        }
        return handle_matrix_action(action, argv[2]);
    }

    // definite integral: integrate "<f>" a b
    if (strcmp(action, "integrate") == 0 && argc == 5) return handle_definite(argv[2], argv[3], argv[4]);

    // polynomial actions require exactly 3 args
    if (argc != 3) {
        print_usage(argv[0]);
        return 1;
    }

    const char* input = argv[2];
    Expr* raw = parse_input(input);
    if (raw && (expr_has_exp(raw) || expr_has_ln(raw) || expr_has_trig(raw))) return handle_exp(action, raw);
    Expr* expr = raw ? simplify(raw) : NULL;
    free_expr(raw);
    if (expr == NULL) {
        fprintf(stderr, "Error: not a polynomial or rational function "
                        "(check syntax; exponents must be integers; division by zero is undefined).\n");
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
            fprintf(stderr, "Error: could not integrate expression: %s.\n", integrate_error());
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
