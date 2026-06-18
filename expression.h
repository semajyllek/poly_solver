#ifndef EXPRESSION_H
#define EXPRESSION_H

#include "polynomial.h"

typedef enum {
    EXPR_CONST,
    EXPR_VAR,
    EXPR_BINOP,
    EXPR_NEG,
    EXPR_POLY,
    EXPR_RATIONAL_FN,
    EXPR_LN
} ExprType;

typedef enum {
    OP_ADD,
    OP_SUB,
    OP_MUL,
    OP_DIV,
    OP_POW
} BinOp;

typedef struct Expr {
    ExprType type;
    union {
        Rational constant;         // EXPR_CONST
        char* var_name;            // EXPR_VAR
        struct {                   // EXPR_BINOP
            BinOp op;
            struct Expr* left;
            struct Expr* right;
        } binop;
        struct Expr* operand;      // EXPR_NEG
        Polynomial* poly;          // EXPR_POLY
        struct {                   // EXPR_RATIONAL_FN
            Polynomial* num;
            Polynomial* den;
        } ratfn;
        struct Expr* ln_arg;       // EXPR_LN: ln|arg|
    };
} Expr;

// constructors
Expr* create_const_rat(Rational value);
Expr* create_const(double value);
Expr* create_var(const char* name);
Expr* create_add(Expr* left, Expr* right);
Expr* create_sub(Expr* left, Expr* right);
Expr* create_mul(Expr* left, Expr* right);
Expr* create_div(Expr* numerator, Expr* denominator);
Expr* create_pow(Expr* base, Expr* exponent);
Expr* create_neg(Expr* operand);
Expr* create_poly_expr(Polynomial* poly);
Expr* create_rational_fn(Polynomial* num, Polynomial* den);
Expr* create_ln(Expr* arg);

// operations
Expr* canonicalize(Expr* expr);
Expr* simplify(Expr* expr);
Expr* factor(Expr* expr);
Expr* differentiate(Expr* expr, const char* variable);
Expr* integrate(Expr* expr, const char* variable);

// output and memory
char* expr_to_string(const Expr* expr);
void print_expr(const Expr* expr);
void free_expr(Expr* expr);

#endif
