#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"          // Include the header for parser functions
#include "expression.h"      // Include the header for expression manipulation
#include "polynomial.h"      // Include the header for polynomial operations
#include "token.h"           // Include token definitions

// Function to display usage instructions
void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s <action> \"<expression>\"\n", program_name);
    fprintf(stderr, "Actions:\n");
    fprintf(stderr, "  integrate   - Integrate the expression\n");
    fprintf(stderr, "  differentiate - Differentiate the expression\n");
    fprintf(stderr, "  factor      - Factor the expression\n");
    fprintf(stderr, "  simplify    - Simplify the expression\n");
}

// Main function
int main(int argc, char* argv[]) {
    // Check for the correct number of arguments
    if (argc != 3) {
        print_usage(argv[0]);
        return 1;
    }

    const char* action = argv[1]; // Action to perform (e.g., integrate, differentiate, factor)
    const char* input = argv[2];   // Mathematical expression input

    // Parse the expression and possibly simplify it
    Expr* expr = parse_and_simplify(input);
    if (expr == NULL) {
        printf("Error parsing expression. Please enter a valid expression.\n");
        return 1;
    }

    // Perform the requested action based on the input
    if (strcmp(action, "factor") == 0) {
        Expr* factored_expr = factor(expr);
        printf("Factored Expression: ");
        print_expr(factored_expr); // Print the factored expression
        printf("\n");
        free_expr(factored_expr);
    } else if (strcmp(action, "differentiate") == 0) {
        Expr* derivative = differentiate(expr, "x");
        printf("Differentiated Expression: ");
        print_expr(derivative); // Print the result of differentiation
        printf("\n");
        free_expr(derivative);
    } else if (strcmp(action, "integrate") == 0) {
        Expr* integral = integrate(expr, "x");
        printf("Integrated Expression: ");
        print_expr(integral); // Print the result of integration
        printf("\n");
        free_expr(integral);
    } else if (strcmp(action, "simplify") == 0) {
        Expr* simplified_expr = simplify(expr);
        printf("Simplified Expression: ");
        print_expr(simplified_expr); // Print the simplified expression
        printf("\n");
        free_expr(simplified_expr);
    } else {
        fprintf(stderr, "Error: Unknown action '%s'.\n", action);
        free_expr(expr);
        return 1;
    }

    // Free the original expression
    free_expr(expr);
    
    return 0;
}
