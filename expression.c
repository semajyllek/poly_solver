#include "expression.h"
#include "polynomial.h"
#include <string.h>


Expr* create_rational(Polynomial* numerator, Polynomial* denominator) {
	RationalExpr* rational = malloc(sizeof(RationalExpr));
	rational->numerator = numerator;
	rational->denominator = denominator;
	return rational;
}

void free_rational(RationalExpr* rational) {
	free_polynomial(rational->numerator);
	free_polynomial(rational->denominator);
	free(rational);
}

void print_rational(const RationalExpr* rational) {
	printf("(");
	print_polynomial(rational->numerator);
	printf(") / (");
	print_polynomial(rational->denominator);
	printf(")");
}

Expr* create_var(const char* name) {
    Expr* expr = malloc(sizeof(Expr));
    expr->type = VAR;
    expr->var_name = strdup(name);
    return expr;
}

Expr* create_const(double value) {
    Expr* expr = malloc(sizeof(Expr));
    expr->type = CONST;
    expr->constant = value;
    return expr;
}

Expr* create_add(Expr* left, Expr* right) {
    Expr* expr = malloc(sizeof(Expr));
    expr->type = ADD;
    expr->add.left = left;
    expr->add.right = right;
    return expr;
}

Expr* create_mul(Expr* left, Expr* right) {
    Expr* expr = malloc(sizeof(Expr));
    expr->type = MUL;
    expr->mul.left = left;
    expr->mul.right = right;
    return expr;
}

Expr* create_pow(Expr* base, Expr* exponent) {
    Expr* expr = malloc(sizeof(Expr));
    expr->type = POW;
    expr->pow.base = base;
    expr->pow.exponent = exponent;
    return expr;
}



Polynomial* polynomial_from_expr(Expr* expr) {
    if (expr->type == POLY) {
        // If the expression is already type POLY
        return expr->polynomial; // Return the associated polynomial
    }

    // Logic to convert other Expr types to Polynomial
    Polynomial* poly = create_polynomial();

    switch (expr->type) {
        case ADD:
            // Assuming both sides of addition are polynomials
            {
                Polynomial* left_poly = polynomial_from_expr(expr->add.left);
                Polynomial* right_poly = polynomial_from_expr(expr->add.right);
                return poly_add(left_poly, right_poly); // Add the two polynomials
            }
        case MUL:
            // Similarly handle multiplication
            {
                Polynomial* left_poly = polynomial_from_expr(expr->mul.left);
                Polynomial* right_poly = polynomial_from_expr(expr->mul.right);
                return poly_multiply(left_poly, right_poly); // Multiply the polynomials
            }
        // Handle other cases as needed...
        default:
            // Return or handle error if it's not convertible
            return NULL;
    }

    return poly; // Return the constructed polynomial if created
}



Expr* simplify(Expr* expr) {
    if (!expr) return NULL;

    switch (expr->type) {
        case VAR:
            return create_var(expr->var_name); // A variable is already simplified

        case CONST:
            return create_const(expr->constant); // A constant is already simplified

        case ADD: {
            Expr* left = simplify(expr->add.left);
            Expr* right = simplify(expr->add.right);
            return create_add(left, right);
        }

        case MUL: {
            Expr* left = simplify(expr->mul.left);
            Expr* right = simplify(expr->mul.right);
            return create_mul(left, right);
        }

        case POW: {
            // Handle power simplification, here we always assume it's in terms of x
            if (strcmp(expr->pow.base->var_name, "x") == 0) { // Default to check against "x"
                double exponent = expr->pow.exponent->constant;
                return create_mul(create_const(1 / (exponent + 1)), create_pow(expr->pow.base, create_const(exponent + 1)));
            }
            return expr; // Just return as is if base is not the variable we want to check
        }

        case RATIONAL: {
            // Simplify the rational expressions
            Polynomial* numerator_poly = expr->rational->numerator;
            Polynomial* denominator_poly = expr->rational->denominator;

            // Calculate the GCD of numerator and denominator
            Polynomial* gcd_poly = poly_gcd(numerator_poly, denominator_poly);

            // Divide both the numerator and denominator by the GCD
            Polynomial* simplified_numerator = poly_divide(numerator_poly, gcd_poly);
            Polynomial* simplified_denominator = poly_divide(denominator_poly, gcd_poly);

            // Create the simplified rational expression
            Expr* simplified = create_rational(simplified_numerator, simplified_denominator);
            free_polynomial(gcd_poly); // Free the GCD polynomial to prevent memory leaks
            return simplified;
        }

        default:
            return expr; // Return the unchanged expression for unrecognized types
    }
}


Expr* integrate(Expr* expr, const char* variable) {
    if (!expr) return NULL; // Check for null expression

    switch (expr->type) {
        case VAR:
            return create_mul(create_const(1.0 / 2), create_pow(expr, create_const(2))); // Integral of x

        case CONST:
            return create_mul(expr, create_var(variable)); // ∫C dx = C*x

        case ADD: {
            Expr* left_integral = integrate(expr->add.left, variable);
            Expr* right_integral = integrate(expr->add.right, variable);
            return create_add(left_integral, right_integral); // Combine integrals
        }

        case MUL: {
            // Here we could implement a more specific treatment for multiplication
            // But for now, handle simple scenarios or leave to user-defined context
            return NULL; 
        }

        case POW: {
            if (strcmp(expr->pow.base->var_name, variable) == 0) {
                double exponent = expr->pow.exponent->constant;
                return create_mul(create_const(1 / (exponent + 1)), create_pow(expr->pow.base, create_const(exponent + 1)));
            }
            return NULL; // If base is not the variable for integration
        }

        case RATIONAL: {
            // For rational expressions, integrate with a dedicated method
            return integrate_rational((RationalExpr*)expr, variable);
        }

        // This assumes you may add support for using poly_integrate when the entire expression is polynomial
        default:
            if (expr->type == RATIONAL || expr->type == POLY) {
                return integrate_poly(expr->polynomial); // Use poly_integrate if this Expr is essentially a polynomial context; not traditionally wrapped
            }
            return NULL; // Handle unrecognized expression types
    }
}
