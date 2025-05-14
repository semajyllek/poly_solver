#include "parser.h"
#include "token.h"
#include "expression.h"
#include "polynomial.h" // For polynomial handling
#include "token_list.c" // For the list data structure in shunting yard algorithm
#include "token_stack.c" // For the stack data structure in shunting yard algorithm


#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Helper functions for parsing
Token* get_next_token(const char** input);


// Define operator precedence
int precedence(TokenType type) {
    switch (type) {
        case TOKEN_PLUS:
        case TOKEN_MINUS:
            return 1;
        case TOKEN_MULTIPLY:
        case TOKEN_DIVIDE:
            return 2;
        case TOKEN_POWER:
            return 3;
        default:
            return 0;
    }
}


// Function to parse individual terms
Expr* parse_term(const char** input) {
    Token* token = get_next_token(input);
    if (token == NULL) return NULL;
    
    Expr* left = NULL;

    if (token->type == TOKEN_NUMBER) {
        left = create_const(token->number);
    } else if (token->type == TOKEN_VARIABLE) {
        left = create_var(token->variable);
        free(token->variable); // Free variable name
    } else if (token->type == TOKEN_LEFT_PAREN) {
        left = parse_expression(input);
        token = get_next_token(input); // Expect matching right parenthesis
        if (token->type != TOKEN_RIGHT_PAREN) {
            free(token);
            free_expr(left);
            return NULL; // Error: unmatched parentheses
        }
    } else {
        free(token);
        return NULL; // Error: unrecognized token
    }

    free(token); // Free the token after processing

    // Check for multiplication or division
    token = get_next_token(input);
    if (token) {
        if (token->type == TOKEN_MULTIPLY) {
            Expr* right = parse_term(input);
            left = create_mul(left, right);
        } else if (token->type == TOKEN_DIVIDE) {
            Expr* right = parse_term(input);
            left = create_div(left, right);
        } else {
            // If not multiplication/division, push back the token
            (*input)--; // Move back input pointer for unprocessed token
        }
    }

    // If the left polynomial was created from an addition operation
    if (left->type == ADD || left->type == MUL) {
        left->type = POLY; // Assign POLY type if it represents polynomial operations
    }

    return left; // Return the constructed term or expression
}

// Modifying parse_expression
Expr* parse_expression(const char** input) {
    Expr* left = parse_term(input);
    if (left == NULL) return NULL; // If there's an error, return null

    Token* token = get_next_token(input);
    
    while (token) {
        if (token->type == TOKEN_PLUS) {
            Expr* right = parse_term(input);
            left = create_add(left, right); // Combine polynomials
            left->type = POLY; // Ensure it retains polynomial type
        } else if (token->type == TOKEN_MINUS) {
            Expr* right = parse_term(input);
            left = create_add(left, create_mul(create_const(-1), right)); // Implement subtraction
        } else {
            // Not an addition or subtraction, break out
            break;
        }
        token = get_next_token(input); // Move to the next token
    }

    return left; // Return the fully constructed expression
}




Expr* build_expression_from_postfix(List* output) {
    Stack stack;
    init_stack(&stack);
    
    ListNode* current = output->head;
    while (current != NULL) {
        Token* token = current->token; // Get current token

        if (token->type == TOKEN_NUMBER) {
            // Create constant expression and push onto stack
            Expr* expr = create_const(token->number);
            push(&stack, expr);
        } else if (token->type == TOKEN_VARIABLE) {
            // Create variable expression and push onto stack
            Expr* expr = create_var(token->variable);
            push(&stack, expr);
        } else {
            // It's an operator, pop operands from the stack
            Expr* right = pop(&stack); // Right operand
            Expr* left = pop(&stack);  // Left operand
            
            Expr* expr = NULL; // Initialize expression
            switch (token->type) {
                case TOKEN_PLUS:
                    expr = create_add(left, right);
                    break;
                case TOKEN_MINUS:
                    expr = create_add(left, create_mul(create_const(-1), right)); // Handle subtraction
                    break;
                case TOKEN_MULTIPLY:
                    expr = create_mul(left, right);
                    break;
                case TOKEN_DIVIDE:
                    expr = create_div(left, right); // Handle division
                    break;
                case TOKEN_POWER:
                    expr = create_pow(left, right); // Handle power
                    break;
            }
            // Optionally check if the resulting expression is polynomial
            if ((left->type == ADD || left->type == MUL) && (right->type == ADD || right->type == MUL)) {
                expr->type = POLY; // Mark the resulting expression as a polynomial
            }
            push(&stack, expr); // Push back the combined expression to stack
        }

        current = current->next; // Iterate over the next token
    }

    // The final expression will be the only remaining element in the stack
    return pop(&stack); 
}


// shunting yard version
Expr* parse_expression_sy(const char** input) {
    Stack operatorStack; // Stack to hold operators
    List output;         // List to hold output in postfix notation
    init_stack(&operatorStack);
    initialize_list(&output);

    Token* token;
    
    while ((token = get_next_token(input)) != NULL) {
        if (token->type == TOKEN_NUMBER) {
            // Push number tokens directly to output
            Expr* numExpr = create_const(token->number);
            add_to_list(&output, numExpr);
        } else if (token->type == TOKEN_VARIABLE) {
            // Push variable tokens directly to output
            Expr* varExpr = create_var(token->variable);
            add_to_list(&output, varExpr);
        } else if (token->type == TOKEN_LEFT_PAREN) {
            // Push left parentheses onto the stack
            push(&operatorStack, token);
        } else if (token->type == TOKEN_RIGHT_PAREN) {
            // Pop from the operator stack to output until a left parenthesis is found
            while (!is_empty(&operatorStack) && peek(&operatorStack)->type != TOKEN_LEFT_PAREN) {
                add_to_list(&output, pop(&operatorStack));
            }
            pop(&operatorStack); // Remove the left parenthesis
        } else if (token->type == TOKEN_PLUS || token->type == TOKEN_MINUS ||
                   token->type == TOKEN_MULTIPLY || token->type == TOKEN_DIVIDE || 
                   token->type == TOKEN_POWER) {
            // Handle operators
            while (!is_empty(&operatorStack) && 
                   precedence(peek(&operatorStack)->type) >= precedence(token->type)) {
                add_to_list(&output, pop(&operatorStack)); // Pop to output before pushing the new operator
            }
            push(&operatorStack, token); // Push current operator onto stack
        }
        free(token); // Free the token memory
    }

    // Pop all remaining operators from the stack to the output
    while (!is_empty(&operatorStack)) {
        add_to_list(&output, pop(&operatorStack));
    }

    // Now, we need to build the expression tree from the output list
    return build_expression_from_postfix(&output);  // Function to convert postfix to expression tree
}




// Function to parse rational expressions
RationalExpr* parse_rational_expression(const char** input) {
    Expr* numerator = parse_expression(input); // Parse the numerator
    Token* token = get_next_token(input); // Expecting a division operator

    if (!token || token->type != TOKEN_DIVIDE) {
        free_expr(numerator);
        return NULL; // Invalid rational expression
    }

    Expr* denominator = parse_expression(input); // Parse the denominator
    if (!denominator) {
        free_expr(numerator);
        return NULL; // Invalid rational expression
    }

    // Create and return the rational expression
    Polynomial* num_poly = polynomial_from_expr(numerator); // Function to convert Expr to Polynomial
    Polynomial* denom_poly = polynomial_from_expr(denominator); // Convert Expr to Polynomial

    return create_rational(num_poly, denom_poly); // Return rational expression
}

// A wrapper to parse and simplify an expression
Expr* parse_and_simplify(const char* input) {
    const char* p = input; // Create a pointer for parsing
    Expr* expression = parse_expression(&p); // Parse the input expression
    if (expression == NULL) {
        printf("Error parsing expression.\n");
        return NULL;
    }

    // Optionally, simplify the resulting expression
    return simplify(expression); 
}

