#include "parser.h"
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
    if (is_unary) return 4;
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

// push operator, first popping higher-precedence ops
static void push_operator(ExprStack* exprs, OpStack* ops, OpEntry entry) {
    int p = precedence(entry.type, entry.is_unary);
    bool ra = is_right_assoc(entry.type, entry.is_unary);

    while (!ostack_empty(ops)) {
        OpEntry top = ostack_peek(ops);
        if (top.type == TOKEN_LEFT_PAREN) break;
        int tp = precedence(top.type, top.is_unary);
        if (tp > p || (tp == p && !ra)) {
            ostack_pop(ops);
            apply_op(exprs, top);
        } else {
            break;
        }
    }
    ostack_push(ops, entry);
}

// true if token type can start an operand (triggers implicit mul detection)
static bool is_operand_start(TokenType t) {
    return t == TOKEN_NUMBER || t == TOKEN_VARIABLE || t == TOKEN_LEFT_PAREN;
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
            OpEntry mul = { TOKEN_MULTIPLY, false };
            push_operator(&exprs, &ops, mul);
        }

        if (tt == TOKEN_NUMBER) {
            estack_push(&exprs, create_const(tok->number));
            prev_was_operand = true;

        } else if (tt == TOKEN_VARIABLE) {
            estack_push(&exprs, create_var(tok->variable));
            free(tok->variable);
            prev_was_operand = true;

        } else if (tt == TOKEN_LEFT_PAREN) {
            ostack_push(&ops, (OpEntry){ TOKEN_LEFT_PAREN, false });
            prev_was_operand = false;

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
                apply_op(&exprs, top);
            }
            if (!found) { success = false; free(tok); break; }
            prev_was_operand = true;

        } else if (tt == TOKEN_PLUS || tt == TOKEN_MINUS ||
                   tt == TOKEN_MULTIPLY || tt == TOKEN_DIVIDE ||
                   tt == TOKEN_POWER) {
            // unary minus: minus when previous was not an operand
            bool unary = (tt == TOKEN_MINUS && !prev_was_operand);
            OpEntry entry = { tt, unary };
            push_operator(&exprs, &ops, entry);
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
        if (top.type == TOKEN_LEFT_PAREN) { success = false; break; }
        apply_op(&exprs, top);
    }

    Expr* result = NULL;
    if (success && exprs.top == 0) {
        result = exprs.items[0];
        exprs.top = -1; // prevent estack_free from freeing the result
    } else {
        estack_free(&exprs);
    }

    ostack_free(&ops);
    if (exprs.top >= 0) free(exprs.items);
    else free(exprs.items);

    return result;
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
