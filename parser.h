#ifndef PARSER_H
#define PARSER_H

#include "token.h"
#include "expression.h"

Expr* parse_expression(const char** input);
Expr* parse_and_simplify(const char* input);
Expr* parse_input(const char* input); // the parse tree, unsimplified (NULL with a message on error)

#endif
