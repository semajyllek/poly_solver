#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>        // Include the assert library
#include "parser.h"       // Include parser functions
#include "expression.h"   // Include expression manipulation
#include "polynomial.h"   // Include polynomial operations
#include "token.h"        // Include token definitions



/*********************************************************************** 
*
*  parse and simplify tests
*
*************************************************************************/

// Test Function for Parsing and Simplifying Expressions
void test_parse_and_simplify_simple_expression() {
    const char* input = "(x + 2) * (5x^3 - 2x^2 + 7)";
    Expr* result = parse_and_simplify(input);
    const char* expected = "5x^4 + 8x^3 + 11x + 14"; // Simplified expression
    char* output = ""; // Convert result to string
    assert(strcmp(output, expected) == 0);
    free_expr(result);
}

void test_parse_and_simplify_polynomial_expression() {
    const char* input = "5x^3 + 2x^2 + 7 + 3x^3 - 2x^2 + 4";
    Expr* result = parse_and_simplify(input);
    const char* expected = "8x^3 + 11"; // Expect simplified output
    char* output = ""; // Convert result to string
    assert(strcmp(output, expected) == 0);
    free_expr(result);
}




void test_only_variables() {
    const char* input = "x + y + z";
    Expr* result = parse_and_simplify(input);
    const char* expected = "x + y + z"; // Should return as is
    char* output = ""; // Convert result to string
    assert(strcmp(output, expected) == 0);
    free_expr(result);
}

void test_negation() {
    const char* input = "-(x + 2)";
    Expr* result = parse_and_simplify(input); 
    const char* expected = "-x - 2"; // Expected output after simplification
    char* output = ""; // Convert result to string
    assert(strcmp(output, expected) == 0);
    free_expr(result);
}

void test_nested_expression() {
    const char* input = "(x + 2) * (3 - (x - 1)) + (x - 4)^2";
    Expr* result = parse_and_simplify(input); 
    const char* expected = "2x^2 - 7x + 9"; // Expected result
    char* output = ""; // Convert result to string
    assert(strcmp(output, expected) == 0);
    free_expr(result);
}

void test_rational_expression() {
    const char* input = "(2x + 3) / (x - 1)";
    RationalExpr* rational_result = parse_rational_expression(&input);
    // Expected: Output could be the rational itself made readable
    // This will depend on the effective representation of rational output.
    // Assuming you'd get expected outputs based on the specific integration/simplification methods
}




/*********************************************************************** 
*
*  algebra tests
*
*************************************************************************/



// Test Function for Polynomial Addition
void test_polynomial_addition() {
    Polynomial* poly1 = malloc(sizeof(Polynomial));
    poly1->head = NULL;
    add_term(poly1, 3, 2); // 3x^2
    add_term(poly1, 4, 1); // 4x
    add_term(poly1, 2, 0); // + 2

    Polynomial* poly2 = malloc(sizeof(Polynomial));
    poly2->head = NULL;
    add_term(poly2, 5, 2); // 5x^2
    add_term(poly2, 3, 1); // 3x
    add_term(poly2, 1, 0); // + 1

    Polynomial* sum = poly_add(poly1, poly2);
    const char* expected_sum = "8x^2 + 7x + 3"; // Expected simplified output
    char* output = ""; // Convert result to string
    assert(strcmp(output, expected_sum) == 0); // Check correctness

    // Free polynomials
    free_polynomial(poly1);
    free_polynomial(poly2);
    free_polynomial(sum);
}

// Test Function for Polynomial Subtraction
void test_polynomial_subtraction() {
    Polynomial* poly1 = malloc(sizeof(Polynomial));
    poly1->head = NULL;
    add_term(poly1, 3, 2); // 3x^2
    add_term(poly1, 4, 1); // 4x
    add_term(poly1, 2, 0); // + 2

    Polynomial* poly2 = malloc(sizeof(Polynomial));
    poly2->head = NULL;
    add_term(poly2, 5, 2); // 5x^2
    add_term(poly2, 3, 1); // 3x
    add_term(poly2, 1, 0); // + 1

    Polynomial* difference = poly_subtract(poly1, poly2);
    const char* expected_diff = "-2x^2 + x + 1"; // Expected simplified output
    char* output = ""; // Convert result to string
    assert(strcmp(output, expected_diff) == 0); // Check correctness

    // Free polynomials
    free_polynomial(poly1);
    free_polynomial(poly2);
    free_polynomial(difference);
}

// Function to test polynomial multiplication
void test_polynomial_multiplication() {
    // Create first polynomial: 3x^2 + 4x + 2
    Polynomial* poly1 = malloc(sizeof(Polynomial));
    poly1->head = NULL;
    add_term(poly1, 3, 2); // 3x^2
    add_term(poly1, 4, 1); // + 4x
    add_term(poly1, 2, 0); // + 2

    // Create second polynomial: 5x^2 + 3x + 1
    Polynomial* poly2 = malloc(sizeof(Polynomial));
    poly2->head = NULL;
    add_term(poly2, 5, 2); // 5x^2
    add_term(poly2, 3, 1); // + 3x
    add_term(poly2, 1, 0); // + 1

    // Perform multiplication
    Polynomial* product = poly_multiply(poly1, poly2);
    const char* expected_product_output = "15x^4 + 22x^3 + 13x^2 + 10x + 2"; // Expected simplified output
    char* product_output = ""; // Convert product to string representation
    assert(strcmp(product_output, expected_product_output) == 0); // Assert if multiplication is correct

    // Free allocated memory
    free_polynomial(poly1);
    free_polynomial(poly2);
    free_polynomial(product);
}

// Function to test polynomial division
void test_polynomial_division() {
    // Create dividend: 6x^3 + 11x^2 + 2
    Polynomial* dividend = malloc(sizeof(Polynomial));
    dividend->head = NULL;
    add_term(dividend, 6, 3); // 6x^3
    add_term(dividend, 11, 2); // + 11x^2
    add_term(dividend, 2, 0); // + 2

    // Create divisor: 3x + 1
    Polynomial* divisor = malloc(sizeof(Polynomial));
    divisor->head = NULL;
    add_term(divisor, 3, 1); // 3x
    add_term(divisor, 1, 0);  // + 1

    // Perform division
    Polynomial* quotient = poly_divide(dividend, divisor);
    const char* expected_quotient_output = "2x^2 + 3x + 1"; // Expected output for quotient
    char* quotient_output = ""; // Convert quotient to string representation
    assert(strcmp(quotient_output, expected_quotient_output) == 0); // Assert if division is correct

    // Free allocated memory
    free_polynomial(dividend);
    free_polynomial(divisor);
    free_polynomial(quotient);
}


void test_polynomial_gcd() {
    // Create first polynomial: x^3 - 2x^2 + x
    Polynomial* poly1 = create_polynomial();
    add_term(poly1, 1, 3); // x^3
    add_term(poly1, -2, 2); // -2x^2
    add_term(poly1, 1, 1); // + x

    // Create second polynomial: x^2 - 1
    Polynomial* poly2 = create_polynomial();
    add_term(poly2, 1, 2); // x^2
    add_term(poly2, -1, 0); // -1
    
    // Perform GCD
    Polynomial* gcd = poly_gcd(poly1, poly2);
    printf("GCD: ");
    print_polynomial(gcd); // Expected output should be x^2 - 1 if correct logic is applied
    printf("\n");

    // Free allocated memory
    free_polynomial(poly1);
    free_polynomial(poly2);
    free_polynomial(gcd);
}






/*********************************************************************** 
*
*  differentiation and integration Tests
*
*************************************************************************/

// Test Function for Differentiation
void test_differentiation() {
    const char* poly_expr = "3x^2 + 4x + 2";
    Expr* expr = parse_and_simplify(poly_expr);
    Expr* derivative = differentiate(expr, "x");
    
    const char* expected_derivative = "6x + 4"; // Expected output for differentiation
    char* output = ""; // Convert result to string
    assert(strcmp(output, expected_derivative) == 0); // Compare outputs
    
    // Free memory
    free_expr(expr);
    free_expr(derivative);
}

// Test Function for Integration
void test_integration() {
    const char* poly_expr = "2x + 3";
    Expr* expr = parse_and_simplify(poly_expr);
    Expr* integral = integrate(expr, "x");
    
    const char* expected_integral_output = "x^2 + 3x + C"; // Expected output for integration (with constant C)
    char* output = ""; // Convert integral result to string representation
    assert(strcmp(output, expected_integral_output) == 0); // Compare outputs
    
    // Free memory
    free_expr(expr);
    free_expr(integral);
}


// This function runs all tests
void run_tests() {
    printf("Running Parsing and Simplification Tests...\n");
    test_parse_and_simplify_simple_expression();
    test_parse_and_simplify_polynomial_expression();
	test_only_variables();
    test_negation();
    test_nested_expression();
	test_rational_expression();

	printf("Running Algebra Tests...\n");
	test_polynomial_addition();
	test_polynomial_subtraction();
	test_polynomial_multiplication();
	test_polynomial_division();
    test_polynomial_gcd();

	printf("Running Differentiation and Integration Tests...\n");
	test_differentiation();
	test_integration();
	
	printf("All tests passed successfully!\n");
	
}

int main() {
    run_tests(); // Execute all tests
    return 0;
}
