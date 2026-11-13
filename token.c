#define _POSIX_C_SOURCE 200809L
#include "token.h"
#include "rational_num.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Function to get the next token from the input string
Token* get_next_token(const char** input) {
    // Skip whitespace
    while (isspace(**input)) (*input)++;

    // Handle end of input
    if (**input == '\0') {
        return NULL; // End of input
    }

    Token* token = (Token*)malloc(sizeof(Token));

    // Handle numbers
    if (isdigit(**input)) {
        token->type = TOKEN_NUMBER;
        // find the literal's extent with the exact parser, then keep its text
        Rational scratch;
        rat_init(&scratch);
        const char* end = rat_parse_literal(&scratch, *input);
        rat_clear(&scratch);
        token->number_text = strndup(*input, (size_t)(end - *input));
        *input = end;
        return token;
    }

    // Handle variables
    if (isalpha(**input)) {
        token->type = TOKEN_VARIABLE;
        token->variable = strndup(*input, 1); // Store the variable character
        (*input)++;
        return token;
    }

    // Handle operators and parentheses
    switch (**input) {
        case '+':
            token->type = TOKEN_PLUS;
            (*input)++;
            break;
        case '-':
            token->type = TOKEN_MINUS;
            (*input)++;
            break;
        case '*':
            token->type = TOKEN_MULTIPLY;
            (*input)++;
            break;
        case '/':
            token->type = TOKEN_DIVIDE;
            (*input)++;
            break;
        case '^':
            token->type = TOKEN_POWER;
            (*input)++;
            break;
        case '(':
            token->type = TOKEN_LEFT_PAREN;
            (*input)++;
            break;
        case ')':
            token->type = TOKEN_RIGHT_PAREN;
            (*input)++;
            break;
        default:
            free(token);
            return NULL; // Unknown character
    }

    return token;
}
