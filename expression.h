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
    EXPR_LN,
    EXPR_ATAN,
    EXPR_SQRT
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
        struct Expr* fn_arg;       // EXPR_ATAN: atan(arg), EXPR_SQRT: sqrt(arg)
    };
} Expr;

// constructors
Expr* create_const_rat(const Rational* value);
Expr* create_const(long value);              // exact integer constant
Expr* create_const_str(const char* literal); // exact decimal literal ("1.25", "2e30")
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
Expr* create_atan(Expr* arg);
Expr* create_sqrt(Expr* arg);

// operations
Expr* canonicalize(Expr* expr);
Expr* simplify(Expr* expr);
Expr* factor(Expr* expr);
Expr* differentiate(Expr* expr, const char* variable);
Expr* integrate(Expr* expr, const char* variable); // integrate.c
const char* integrate_error(void);                  // reason for the last NULL from integrate

// numeric value at x (for verification); *ok becomes false on a pole or domain error
double expr_eval(const Expr* expr, double x, bool* ok);

// output and memory
char* expr_to_string(const Expr* expr);
void print_expr(const Expr* expr);
void free_expr(Expr* expr);

#endif
