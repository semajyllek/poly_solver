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

    // functions: exp( ln( sin( cos( tan( sqrt( (only when the name is followed by an opening parenthesis)
    static const struct { const char* name; FunctionName fn; } FUNCTIONS[] = {
        { "exp", FN_EXP }, { "ln", FN_LN }, { "sin", FN_SIN }, { "cos", FN_COS }, { "tan", FN_TAN }, { "sqrt", FN_SQRT } };
    for (size_t i = 0; i < sizeof FUNCTIONS / sizeof FUNCTIONS[0]; i++) {
        size_t n = strlen(FUNCTIONS[i].name);
        if (strncmp(*input, FUNCTIONS[i].name, n) != 0) continue;   // also stops at the string's end
        const char* after = *input + n;
        while (isspace((unsigned char)*after)) after++;
        if (*after == '(') {
            token->type = TOKEN_FUNCTION;
            token->function = FUNCTIONS[i].fn;
            *input = after;                           // the '(' is the next token
            return token;
        }
    }

    // the variable: x only (any other letter is an error, not a silent second name for x)
    if (**input == 'x') {
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
