#ifndef PARSER_H
#define PARSER_H

#include "token.h"        // Include the token definitions
#include "expression.h"   // Include definitions for expression structures

// Function prototypes for parsing
Expr* parse_expression(const char** input);
Polynomial* parse_polynomial(const char** input); // Function specifically for parsing polynomials
RationalExpr* parse_rational_expression(const char** input); // Function for rational expressions

// Function to parse and simplify an expression
Expr* parse_and_simplify(const char* input);

#endif // PARSER_H
