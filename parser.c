#include "parser.h"
#include <ctype.h>
#include "token.h"
#include "expression.h"
#include "polynomial.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// --- operand stack (Expr*) ---

typedef struct {
    Expr** items;
    int top;
    int cap;
} ExprStack;

static ExprStack estack_new(void) {
    ExprStack s = { malloc(16 * sizeof(Expr*)), -1, 16 };
    return s;
}

static void estack_push(ExprStack* s, Expr* e) {
    if (s->top + 1 >= s->cap) {
        s->cap *= 2;
        s->items = realloc(s->items, s->cap * sizeof(Expr*));
    }
    s->items[++s->top] = e;
}

static Expr* estack_pop(ExprStack* s) {
    if (s->top < 0) return NULL;
    return s->items[s->top--];
}

static void estack_free(ExprStack* s) {
    for (int i = 0; i <= s->top; i++)
        free_expr(s->items[i]);
    free(s->items);
}

// --- operator stack ---

typedef struct {
    TokenType type;
    bool is_unary;
    FunctionName fn; // for TOKEN_FUNCTION
} OpEntry;

typedef struct {
    OpEntry* items;
    int top;
    int cap;
} OpStack;

static OpStack ostack_new(void) {
    OpStack s = { malloc(16 * sizeof(OpEntry)), -1, 16 };
    return s;
}

static void ostack_push(OpStack* s, OpEntry e) {
    if (s->top + 1 >= s->cap) {
        s->cap *= 2;
        s->items = realloc(s->items, s->cap * sizeof(OpEntry));
    }
    s->items[++s->top] = e;
}

static OpEntry ostack_pop(OpStack* s) {
    return s->items[s->top--];
}

static OpEntry ostack_peek(OpStack* s) {
    return s->items[s->top];
}

static bool ostack_empty(OpStack* s) {
    return s->top < 0;
}

static void ostack_free(OpStack* s) {
    free(s->items);
}

// --- precedence and associativity ---

static int precedence(TokenType type, bool is_unary) {
    // unary minus binds like ^ (both right-associative): -x^2 = -(x^2), x^-2 = x^(-2)
    if (is_unary) return 3;
    switch (type) {
        case TOKEN_PLUS: case TOKEN_MINUS: return 1;
        case TOKEN_MULTIPLY: case TOKEN_DIVIDE: return 2;
        case TOKEN_POWER: return 3;
        default: return 0;
    }
}

static bool is_right_assoc(TokenType type, bool is_unary) {
    return is_unary || type == TOKEN_POWER;
}

// apply operator on top of op stack to operands on expr stack
static bool apply_op(ExprStack* exprs, OpEntry op) {
    if (op.is_unary) {
        Expr* operand = estack_pop(exprs);
        if (!operand) return false;
        estack_push(exprs, create_neg(operand));
        return true;
    }

    Expr* right = estack_pop(exprs);
    Expr* left = estack_pop(exprs);
    if (!left || !right) {
        free_expr(left);
        free_expr(right);
        return false;
    }

    Expr* result = NULL;
    switch (op.type) {
        case TOKEN_PLUS:     result = create_add(left, right); break;
        case TOKEN_MINUS:    result = create_sub(left, right); break;
        case TOKEN_MULTIPLY: result = create_mul(left, right); break;
        case TOKEN_DIVIDE:   result = create_div(left, right); break;
        case TOKEN_POWER:    result = create_pow(left, right); break;
        default:
            free_expr(left);
            free_expr(right);
            return false;
    }
    estack_push(exprs, result);
    return true;
}

// push operator, first popping higher-precedence ops; false if an operand was missing
static bool push_operator(ExprStack* exprs, OpStack* ops, OpEntry entry) {
    bool ok = true;
    int p = precedence(entry.type, entry.is_unary);
    bool ra = is_right_assoc(entry.type, entry.is_unary);

    while (!ostack_empty(ops)) {
        OpEntry top = ostack_peek(ops);
        if (top.type == TOKEN_LEFT_PAREN) break;
        int tp = precedence(top.type, top.is_unary);
        if (tp > p || (tp == p && !ra)) {
            ostack_pop(ops);
            ok = apply_op(exprs, top) && ok;
        } else {
            break;
        }
    }
    ostack_push(ops, entry);
    return ok;
}

// true if token type can start an operand (triggers implicit mul detection)
static bool is_operand_start(TokenType t) {
    return t == TOKEN_NUMBER || t == TOKEN_VARIABLE || t == TOKEN_LEFT_PAREN || t == TOKEN_FUNCTION;
}

// apply the function waiting below a just-closed parenthesis, if there is one
static bool apply_function(ExprStack* exprs, OpStack* ops) {
    if (ostack_empty(ops) || ostack_peek(ops).type != TOKEN_FUNCTION) return true;
    OpEntry f = ostack_pop(ops);
    Expr* arg = estack_pop(exprs);
    if (!arg) return false;
    estack_push(exprs, f.fn == FN_EXP ? create_exp(arg) : create_ln(arg));
    return true;
}

// --- shunting-yard parser ---

Expr* parse_expression(const char** input) {
    ExprStack exprs = estack_new();
    OpStack ops = ostack_new();
    bool prev_was_operand = false;
    bool success = true;

    while (1) {
        const char* saved = *input;
        Token* tok = get_next_token(input);
        if (!tok) break;

        TokenType tt = tok->type;

        // implicit multiplication: operand followed by operand-start
        if (prev_was_operand && is_operand_start(tt)) {
            OpEntry mul = { TOKEN_MULTIPLY, false, FN_EXP };
            success = push_operator(&exprs, &ops, mul) && success;
        }

        if (tt == TOKEN_NUMBER) {
            estack_push(&exprs, create_const_str(tok->number_text));
            free(tok->number_text);
            prev_was_operand = true;

        } else if (tt == TOKEN_VARIABLE) {
            estack_push(&exprs, create_var(tok->variable));
            free(tok->variable);
            prev_was_operand = true;

        } else if (tt == TOKEN_LEFT_PAREN) {
            ostack_push(&ops, (OpEntry){ TOKEN_LEFT_PAREN, false, FN_EXP });
            prev_was_operand = false;

        } else if (tt == TOKEN_FUNCTION) {
            ostack_push(&ops, (OpEntry){ TOKEN_FUNCTION, false, tok->function });
            prev_was_operand = false;                   // its '(' follows (the lexer checked)

        } else if (tt == TOKEN_RIGHT_PAREN) {
            // pop until matching left paren
            bool found = false;
            while (!ostack_empty(&ops)) {
                OpEntry top = ostack_peek(&ops);
                if (top.type == TOKEN_LEFT_PAREN) {
                    ostack_pop(&ops);
                    found = true;
                    break;
                }
                ostack_pop(&ops);
                success = apply_op(&exprs, top) && success;
            }
            if (!found) { success = false; free(tok); break; }
            success = apply_function(&exprs, &ops) && success;
            prev_was_operand = true;

        } else if (tt == TOKEN_PLUS || tt == TOKEN_MINUS ||
                   tt == TOKEN_MULTIPLY || tt == TOKEN_DIVIDE ||
                   tt == TOKEN_POWER) {
            // unary minus: minus when previous was not an operand
            bool unary = (tt == TOKEN_MINUS && !prev_was_operand);
            OpEntry entry = { tt, unary, FN_EXP };
            success = push_operator(&exprs, &ops, entry) && success;
            prev_was_operand = false;

        } else {
            // unknown token — put back and stop
            *input = saved;
            free(tok);
            break;
        }

        free(tok);
    }

    // pop remaining operators
    while (!ostack_empty(&ops)) {
        OpEntry top = ostack_pop(&ops);
        if (top.type == TOKEN_LEFT_PAREN || top.type == TOKEN_FUNCTION) { success = false; break; }
        success = apply_op(&exprs, top) && success;
    }

    Expr* result = NULL;
    if (success && exprs.top == 0) {
        result = exprs.items[0];
        free(exprs.items); // the result now belongs to the caller
    } else {
        estack_free(&exprs); // frees remaining expressions and the array, exactly once
    }
    ostack_free(&ops);
    return result;
}

Expr* parse_and_simplify(const char* input) {
    Expr* expression = parse_input(input);
    if (!expression) return NULL;
    Expr* result = simplify(expression);
    free_expr(expression);
    return result;
}

Expr* parse_input(const char* input) {
    const char* p = input;
    Expr* expression = parse_expression(&p);
    while (expression && isspace((unsigned char)*p)) p++;
    if (expression && *p != '\0') { // never silently ignore trailing input such as "= 4"
        fprintf(stderr, "Error: unexpected '%c' at position %ld.\n", *p, (long)(p - input) + 1);
        free_expr(expression);
        return NULL;
    }
    if (expression == NULL) {
        fprintf(stderr, "Error parsing expression.\n");
        return NULL;
    }
    return expression;
}
