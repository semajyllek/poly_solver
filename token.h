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
    TOKEN_END
} TokenType;

typedef struct {
    TokenType type; // Type of token
    union {
        double number;  // For constant numeric values
        char* variable; // For variable names
    };
} Token;

// Tokens related functions
Token* get_next_token(const char** input);

#endif // TOKEN_H
