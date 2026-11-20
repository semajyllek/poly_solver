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
    EXPR_SQRT,
    EXPR_ROOTSUM,
    EXPR_TRIG,
    EXPR_EXP,
    EXPR_SIN,
    EXPR_COS
} ExprType;

// exact trig constants at rational multiples of π (EXPR_TRIG); TRIG_NONE is angle·π itself
typedef enum { TRIG_NONE, TRIG_COS, TRIG_SIN, TRIG_COT, TRIG_CSC } TrigKind;

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
        struct Expr* ln_arg;       // EXPR_LN: ln(abs(arg))
        struct Expr* fn_arg;       // EXPR_ATAN: atan(arg), EXPR_SQRT: sqrt(arg), EXPR_EXP: exp(arg), EXPR_SIN/COS: sin/cos(arg)
        struct {                   // EXPR_ROOTSUM: Σ over roots t of q of t·ln(S(t, x))
            Polynomial* q;         //   univariate in t, irreducible over Q
            Polynomial** s;        //   S = Σ s[k](t)·x^k, each s[k] of degree < deg q (NULL with c)
            int s_deg;
            struct Expr** c;       //   or S = Σ t^i·c[i], c[i] an expression in x (NULL for 0), i < deg q:
        } rootsum;                 //   over exp/ln/tan towers, where S is a polynomial in θ
        struct {                   // EXPR_TRIG: cos/sin/cot/csc(angle·π), or angle·π (TRIG_NONE)
            TrigKind kind;
            Rational angle;
        } trig;
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
Expr* expr_append_sum(Expr* acc, Expr* e);          // acc + e, e's terms spliced in; takes both, acc may be NULL
Expr* create_atan(Expr* arg);
Expr* create_exp(Expr* arg);
Expr* create_sin(Expr* arg);
Expr* create_cos(Expr* arg);
Expr* create_sqrt(Expr* arg);
Expr* create_rootsum(Polynomial* q, Polynomial** s, int s_deg); // takes ownership
Expr* create_rootsum_terms(Polynomial* q, Expr** c);             // S = Σ t^i·c[i]; takes ownership
Expr* create_trig(TrigKind kind, const Rational* angle);        // trig(angle·π)
Expr* create_pi(void);                                          // π
Expr* expr_copy(const Expr* e);                                 // deep copy
bool  expr_has_rootsum(const Expr* e);                          // a RootSum among the terms
bool  expr_has_exp(const Expr* e);                              // an exp(...) anywhere in e
bool  expr_has_ln(const Expr* e);                               // an ln(...) anywhere in e
bool  expr_has_trig(const Expr* e);                             // a sin(...) or cos(...) anywhere in e

// operations
Expr* canonicalize(Expr* expr);
Expr* simplify(Expr* expr);
Expr* factor(Expr* expr);
Expr* differentiate(Expr* expr, const char* variable);
Expr* integrate(Expr* expr, const char* variable); // integrate.c
Expr* integrate_rootsum(Expr* expr, const char* variable); // same, but RootSum is never expanded
const char* integrate_error(void);                  // reason for the last NULL from integrate

// numeric value at x (for verification); *ok becomes false on a pole or domain error
double expr_eval(const Expr* expr, double x, bool* ok);

// output and memory
char* expr_to_string(const Expr* expr);
void print_expr(const Expr* expr);
void free_expr(Expr* expr);

#endif
