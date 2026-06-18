#ifndef PARSER_H
#define PARSER_H

#include "token.h"
#include "expression.h"

Expr* parse_expression(const char** input);
Expr* parse_and_simplify(const char* input);

#endif
