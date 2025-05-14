#ifndef EXPRESSION_H
#define EXPRESSION_H

typedef struct Polynomial Polynomial; // Forward declaration

typedef struct RationalExpr {
    Polynomial* numerator;   // Pointer to the numerator polynomial
    Polynomial* denominator; // Pointer to the denominator polynomial
} RationalExpr;

typedef enum {
    VAR,
    CONST,
    ADD,
    SUBTRACT,
    MUL,
    DIVIDE,
    POW,
    POLY,       // For explicitly indicating polynomial expressions
    RATIONAL
} ExprType;

typedef struct Expr {
    ExprType type; // Type of expression (VAR, CONST, ADD, etc.)
    union {
        char* var_name;      // Variable name for VAR type
        double constant;     // Constant value for CONST type
        struct {            // Addition structure
            struct Expr* left;    // Left side of addition
            struct Expr* right;   // Right side of addition
        } add;
        struct {            // Subtraction structure
            struct Expr* left;    // Left side of subtraction
            struct Expr* right;   // Right side of subtraction
        } subtract;
        struct {            // Multiplication structure
            struct Expr* left;    // Left side of multiplication
            struct Expr* right;   // Right side of multiplication
        } mul;
        struct {            // Division structure
            struct Expr* numerator;   // Numerator for division
            struct Expr* denominator; // Denominator for division
        } div;
        struct {            // Power structure
            struct Expr* base;     // Base for power operation
            struct Expr* exponent; // Exponent for power operation
        } pow;
		Polynomial* polynomial; // Pointer for POLY expressions
        RationalExpr* rational; // Pointer for RATIONAL expressions
    };
} Expr;

// Function prototypes
Expr* simplify(Expr* expr);
Expr* differentiate(Expr* expr, const char* variable);
Expr* integrate(Expr* expr, const char* variable);
void free_expr(Expr* expr);
void print_expr(const Expr* expr);

#endif // EXPRESSION_H
