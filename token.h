#ifndef TOKEN_H
#define TOKEN_H

#include <stdbool.h>

typedef enum {
    TOKEN_NUMBER,
    TOKEN_VARIABLE,
    TOKEN_PLUS,
    TOKEN_MINUS,
    TOKEN_MULTIPLY,
    TOKEN_DIVIDE,
    TOKEN_POWER,
    TOKEN_LEFT_PAREN,
    TOKEN_RIGHT_PAREN,
    TOKEN_FUNCTION, // exp( or ln( : the name, with its '(' still to come
    TOKEN_END
} TokenType;

typedef enum { FN_EXP, FN_LN } FunctionName;

typedef struct {
    TokenType type; // Type of token
    union {
        char* number_text; // exact literal text, e.g. "3", "1.25", "2e10"
        char* variable; // For variable names
        FunctionName function;
    };
} Token;

// Tokens related functions
Token* get_next_token(const char** input);

#endif // TOKEN_H
