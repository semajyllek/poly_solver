#include "parser.h"
#include "token.h"
#include "expression.h"
#include "polynomial.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


static Expr* parse_term(const char** input);

static Expr* parse_atom(const char** input) {
    Token* token = get_next_token(input);
    if (token == NULL) return NULL;

    Expr* result = NULL;

    if (token->type == TOKEN_NUMBER) {
        result = create_const(token->number);
    } else if (token->type == TOKEN_VARIABLE) {
        result = create_var(token->variable);
        free(token->variable);
    } else if (token->type == TOKEN_LEFT_PAREN) {
        free(token);
        result = parse_expression(input);
        Token* rparen = get_next_token(input);
        if (!rparen || rparen->type != TOKEN_RIGHT_PAREN) {
            free(rparen);
            free_expr(result);
            return NULL;
        }
        free(rparen);
        return result;
    } else {
        free(token);
        return NULL;
    }

    free(token);
    return result;
}

static Expr* parse_term(const char** input) {
    const char* saved = *input;
    Token* tok = get_next_token(input);

    // handle unary minus
    if (tok && tok->type == TOKEN_MINUS) {
        free(tok);
        Expr* operand = parse_term(input);
        if (!operand) return NULL;
        return create_neg(operand);
    }

    // not unary minus — restore and parse atom
    *input = saved;
    if (tok) free(tok);

    Expr* base = parse_atom(input);
    if (!base) return NULL;

    // check for implicit multiplication: number followed by variable or paren
    saved = *input;
    tok = get_next_token(input);
    if (tok && (tok->type == TOKEN_VARIABLE || tok->type == TOKEN_LEFT_PAREN)) {
        // implicit multiplication
        *input = saved;
        free(tok);
        Expr* right = parse_term(input);
        if (right)
            base = create_mul(base, right);
        return base;
    }

    // check for power
    if (tok && tok->type == TOKEN_POWER) {
        free(tok);
        Expr* exponent = parse_term(input);
        if (!exponent) return base;

        // check for implicit mul after power: x^2x should be x^2 * x
        saved = *input;
        tok = get_next_token(input);
        if (tok && (tok->type == TOKEN_VARIABLE || tok->type == TOKEN_LEFT_PAREN || tok->type == TOKEN_NUMBER)) {
            *input = saved;
            free(tok);
            Expr* pow_expr = create_pow(base, exponent);
            Expr* right = parse_term(input);
            if (right)
                return create_mul(pow_expr, right);
            return pow_expr;
        }
        if (tok) { *input = saved; free(tok); }
        return create_pow(base, exponent);
    }

    // not an operator — restore
    if (tok) { *input = saved; free(tok); }
    return base;
}

Expr* parse_expression(const char** input) {
    Expr* left = parse_term(input);
    if (left == NULL) return NULL;

    while (1) {
        const char* saved = *input;
        Token* token = get_next_token(input);
        if (!token) break;

        if (token->type == TOKEN_PLUS) {
            free(token);
            Expr* right = parse_term(input);
            if (!right) { free_expr(left); return NULL; }
            left = create_add(left, right);
        } else if (token->type == TOKEN_MINUS) {
            free(token);
            Expr* right = parse_term(input);
            if (!right) { free_expr(left); return NULL; }
            left = create_sub(left, right);
        } else if (token->type == TOKEN_MULTIPLY) {
            free(token);
            Expr* right = parse_term(input);
            if (!right) { free_expr(left); return NULL; }
            left = create_mul(left, right);
        } else if (token->type == TOKEN_DIVIDE) {
            free(token);
            Expr* right = parse_term(input);
            if (!right) { free_expr(left); return NULL; }
            left = create_div(left, right);
        } else {
            *input = saved;
            free(token);
            break;
        }
    }

    return left;
}

Expr* parse_and_simplify(const char* input) {
    const char* p = input;
    Expr* expression = parse_expression(&p);
    if (expression == NULL) {
        fprintf(stderr, "Error parsing expression.\n");
        return NULL;
    }
    Expr* result = simplify(expression);
    free_expr(expression);
    return result;
}
