#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "rational_num.h"
#include "parser.h"
#include "expression.h"
#include "polynomial.h"
#include "factor.h"
#include "solver.h"
#include "token.h"


/***********************************************************************
*  rational number tests
*************************************************************************/

void test_rat_normalization() {
    Rational r1 = rat_create(6, 4);
    assert(r1.num == 3 && r1.den == 2);

    Rational r2 = rat_create(-6, -4);
    assert(r2.num == 3 && r2.den == 2);

    Rational r3 = rat_create(6, -4);
    assert(r3.num == -3 && r3.den == 2);

    Rational r4 = rat_create(0, 5);
    assert(r4.num == 0 && r4.den == 1);

    printf("  test_rat_normalization PASSED\n");
}

void test_rat_arithmetic() {
    // 1/2 + 1/3 = 5/6
    Rational sum = rat_add(rat_create(1, 2), rat_create(1, 3));
    assert(sum.num == 5 && sum.den == 6);

    // 3/4 - 1/4 = 1/2
    Rational diff = rat_sub(rat_create(3, 4), rat_create(1, 4));
    assert(diff.num == 1 && diff.den == 2);

    // 2/3 * 3/4 = 1/2
    Rational prod = rat_mul(rat_create(2, 3), rat_create(3, 4));
    assert(prod.num == 1 && prod.den == 2);

    // (2/3) / (4/5) = 10/12 = 5/6
    Rational quot = rat_div(rat_create(2, 3), rat_create(4, 5));
    assert(quot.num == 5 && quot.den == 6);

    // neg
    Rational neg = rat_neg(rat_create(3, 7));
    assert(neg.num == -3 && neg.den == 7);

    printf("  test_rat_arithmetic PASSED\n");
}

void test_rat_pow() {
    // (2/3)^3 = 8/27
    Rational r = rat_pow_int(rat_create(2, 3), 3);
    assert(r.num == 8 && r.den == 27);

    // (2/3)^0 = 1
    Rational r0 = rat_pow_int(rat_create(2, 3), 0);
    assert(r0.num == 1 && r0.den == 1);

    // (2/3)^-1 = 3/2
    Rational rn = rat_pow_int(rat_create(2, 3), -1);
    assert(rn.num == 3 && rn.den == 2);

    printf("  test_rat_pow PASSED\n");
}

void test_rat_comparison() {
    assert(rat_eq(rat_create(2, 4), rat_create(1, 2)));
    assert(!rat_eq(rat_create(1, 2), rat_create(1, 3)));
    assert(rat_is_zero(rat_create(0, 42)));
    assert(rat_is_one(rat_create(5, 5)));
    assert(rat_is_negative(rat_create(-1, 3)));
    assert(!rat_is_negative(rat_create(1, 3)));

    // cmp: 1/2 < 2/3
    assert(rat_cmp(rat_create(1, 2), rat_create(2, 3)) < 0);
    // cmp: 3/4 > 1/2
    assert(rat_cmp(rat_create(3, 4), rat_create(1, 2)) > 0);
    // cmp: 1/2 == 2/4
    assert(rat_cmp(rat_create(1, 2), rat_create(2, 4)) == 0);

    printf("  test_rat_comparison PASSED\n");
}

void test_rat_to_string() {
    char* s1 = rat_to_string(rat_create(3, 1));
    assert(strcmp(s1, "3") == 0);
    free(s1);

    char* s2 = rat_to_string(rat_create(3, 5));
    assert(strcmp(s2, "3/5") == 0);
    free(s2);

    char* s3 = rat_to_string(rat_create(-1, 2));
    assert(strcmp(s3, "-1/2") == 0);
    free(s3);

    printf("  test_rat_to_string PASSED\n");
}

void test_rat_to_double() {
    assert(fabs(rat_to_double(rat_create(1, 2)) - 0.5) < 1e-15);
    assert(fabs(rat_to_double(rat_create(1, 3)) - 0.333333333333333) < 1e-10);

    printf("  test_rat_to_double PASSED\n");
}


/***********************************************************************
*  helper: shorthand for creating Rational from int
*************************************************************************/

#define R(n) rat_from_int(n)

/***********************************************************************
*  polynomial arithmetic tests
*************************************************************************/

void test_polynomial_addition() {
    Polynomial* poly1 = create_polynomial();
    add_term(poly1, R(3), 2);
    add_term(poly1, R(4), 1);
    add_term(poly1, R(2), 0);

    Polynomial* poly2 = create_polynomial();
    add_term(poly2, R(5), 2);
    add_term(poly2, R(3), 1);
    add_term(poly2, R(1), 0);

    Polynomial* sum = poly_add(poly1, poly2);

    // expect 8x^2 + 7x + 3
    Term* t = sum->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(8)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(7)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(3)));
    assert(t->next == NULL);

    free_polynomial(poly1);
    free_polynomial(poly2);
    free_polynomial(sum);
    printf("  test_polynomial_addition PASSED\n");
}

void test_polynomial_subtraction() {
    Polynomial* poly1 = create_polynomial();
    add_term(poly1, R(3), 2);
    add_term(poly1, R(4), 1);
    add_term(poly1, R(2), 0);

    Polynomial* poly2 = create_polynomial();
    add_term(poly2, R(5), 2);
    add_term(poly2, R(3), 1);
    add_term(poly2, R(1), 0);

    Polynomial* diff = poly_subtract(poly1, poly2);

    // expect -2x^2 + x + 1
    Term* t = diff->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(-2)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    free_polynomial(poly1);
    free_polynomial(poly2);
    free_polynomial(diff);
    printf("  test_polynomial_subtraction PASSED\n");
}

void test_polynomial_multiplication() {
    // (3x^2 + 4x + 2) * (5x^2 + 3x + 1)
    Polynomial* poly1 = create_polynomial();
    add_term(poly1, R(3), 2);
    add_term(poly1, R(4), 1);
    add_term(poly1, R(2), 0);

    Polynomial* poly2 = create_polynomial();
    add_term(poly2, R(5), 2);
    add_term(poly2, R(3), 1);
    add_term(poly2, R(1), 0);

    Polynomial* product = poly_multiply(poly1, poly2);

    // expect 15x^4 + 29x^3 + 25x^2 + 10x + 2
    Term* t = product->head;
    assert(t && t->exponent == 4 && rat_eq(t->coeff, R(15)));
    t = t->next;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(29)));
    t = t->next;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(25)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(10)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(2)));
    assert(t->next == NULL);

    free_polynomial(poly1);
    free_polynomial(poly2);
    free_polynomial(product);
    printf("  test_polynomial_multiplication PASSED\n");
}

void test_polynomial_division() {
    // (x^2 - 1) / (x - 1) = x + 1
    Polynomial* dividend = create_polynomial();
    add_term(dividend, R(1), 2);
    add_term(dividend, R(-1), 0);

    Polynomial* divisor = create_polynomial();
    add_term(divisor, R(1), 1);
    add_term(divisor, R(-1), 0);

    Polynomial* quotient = poly_divide(dividend, divisor);

    // expect x + 1
    Term* t = quotient->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    free_polynomial(dividend);
    free_polynomial(divisor);
    free_polynomial(quotient);
    printf("  test_polynomial_division PASSED\n");
}

void test_polynomial_divmod() {
    // (x^3 + 2x^2 + 3x + 4) / (x + 1)
    // quotient = x^2 + x + 2, remainder = 2
    Polynomial* dividend = create_polynomial();
    add_term(dividend, R(1), 3);
    add_term(dividend, R(2), 2);
    add_term(dividend, R(3), 1);
    add_term(dividend, R(4), 0);

    Polynomial* divisor = create_polynomial();
    add_term(divisor, R(1), 1);
    add_term(divisor, R(1), 0);

    PolyDivResult dr = poly_divmod(dividend, divisor);

    Term* t = dr.quotient->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(2)));
    assert(t->next == NULL);

    t = dr.remainder->head;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(2)));
    assert(t->next == NULL);

    free_polynomial(dividend);
    free_polynomial(divisor);
    free_polynomial(dr.quotient);
    free_polynomial(dr.remainder);
    printf("  test_polynomial_divmod PASSED\n");
}

void test_polynomial_gcd() {
    // gcd(x^2 - 1, x^2 - 2x + 1) = x - 1
    // x^2-1 = (x-1)(x+1), x^2-2x+1 = (x-1)^2
    Polynomial* poly1 = create_polynomial();
    add_term(poly1, R(1), 2);
    add_term(poly1, R(-1), 0);

    Polynomial* poly2 = create_polynomial();
    add_term(poly2, R(1), 2);
    add_term(poly2, R(-2), 1);
    add_term(poly2, R(1), 0);

    // verify inputs survive the call (regression test for old destructive bug)
    Polynomial* poly1_copy = poly_copy(poly1);

    Polynomial* gcd = poly_gcd(poly1, poly2);

    // inputs should be unmodified
    assert(poly1->head != NULL);
    Term* t1 = poly1->head;
    assert(t1->exponent == 2 && rat_eq(t1->coeff, R(1)));

    // gcd should be monic (x - 1)
    Term* t = gcd->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(-1)));
    assert(t->next == NULL);

    free_polynomial(poly1);
    free_polynomial(poly1_copy);
    free_polynomial(poly2);
    free_polynomial(gcd);
    printf("  test_polynomial_gcd PASSED\n");
}

void test_polynomial_derivative() {
    // d/dx(3x^2 + 4x + 2) = 6x + 4
    Polynomial* poly = create_polynomial();
    add_term(poly, R(3), 2);
    add_term(poly, R(4), 1);
    add_term(poly, R(2), 0);

    Polynomial* deriv = poly_derivative(poly);

    Term* t = deriv->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(6)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(4)));
    assert(t->next == NULL);

    free_polynomial(poly);
    free_polynomial(deriv);
    printf("  test_polynomial_derivative PASSED\n");
}

void test_polynomial_integral() {
    // integral(2x + 3) = x^2 + 3x
    Polynomial* poly = create_polynomial();
    add_term(poly, R(2), 1);
    add_term(poly, R(3), 0);

    Polynomial* integral = poly_integral(poly);

    Term* t = integral->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(3)));
    assert(t->next == NULL);

    free_polynomial(poly);
    free_polynomial(integral);
    printf("  test_polynomial_integral PASSED\n");
}

void test_polynomial_evaluate() {
    // 3x^2 + 4x + 2 at x=2 => 22
    Polynomial* poly = create_polynomial();
    add_term(poly, R(3), 2);
    add_term(poly, R(4), 1);
    add_term(poly, R(2), 0);

    double result = poly_evaluate(poly, 2.0);
    assert(fabs(result - 22.0) < 1e-10);

    free_polynomial(poly);
    printf("  test_polynomial_evaluate PASSED\n");
}

void test_polynomial_degree() {
    Polynomial* poly = create_polynomial();
    add_term(poly, R(3), 5);
    add_term(poly, R(1), 2);
    add_term(poly, R(7), 0);

    assert(polynomial_degree(poly) == 5);

    free_polynomial(poly);
    printf("  test_polynomial_degree PASSED\n");
}

void test_poly_to_string() {
    Polynomial* poly = create_polynomial();
    add_term(poly, R(3), 2);
    add_term(poly, R(-4), 1);
    add_term(poly, R(5), 0);

    char* s = poly_to_string(poly);
    assert(strcmp(s, "3x^2 - 4x + 5") == 0);
    free(s);

    // rational coefficients
    Polynomial* poly2 = create_polynomial();
    add_term(poly2, rat_create(1, 2), 2);
    add_term(poly2, rat_create(-3, 4), 0);

    char* s2 = poly_to_string(poly2);
    assert(strcmp(s2, "1/2x^2 - 3/4") == 0);
    free(s2);

    free_polynomial(poly);
    free_polynomial(poly2);
    printf("  test_poly_to_string PASSED\n");
}

void test_poly_copy() {
    Polynomial* poly = create_polynomial();
    add_term(poly, R(3), 2);
    add_term(poly, R(1), 0);

    Polynomial* copy = poly_copy(poly);

    // modify original, copy should be unaffected
    add_term(poly, R(5), 1);

    Term* t = copy->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(3)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    free_polynomial(poly);
    free_polynomial(copy);
    printf("  test_poly_copy PASSED\n");
}

/***********************************************************************
*  parser tests
*************************************************************************/

void test_parse_implicit_mul() {
    // "5x^3" should parse and simplify to 5x^3
    Expr* e = parse_and_simplify("5x^3");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(5)));
    assert(t->next == NULL);
    free_expr(e);
    printf("  test_parse_implicit_mul PASSED\n");
}

void test_parse_unary_minus() {
    // "-x" should parse to -1x
    Expr* e = parse_and_simplify("-x");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(-1)));
    assert(t->next == NULL);
    free_expr(e);
    printf("  test_parse_unary_minus PASSED\n");
}

void test_parse_and_simplify_polynomial() {
    // "5x^3 + 3x^3 - 2x + 4 + 7" => 8x^3 - 2x + 11
    Expr* e = parse_and_simplify("5x^3 + 3x^3 - 2x + 4 + 7");
    assert(e && e->type == EXPR_POLY);

    Term* t = e->poly->head;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(8)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(-2)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(11)));
    assert(t->next == NULL);

    free_expr(e);
    printf("  test_parse_and_simplify_polynomial PASSED\n");
}

/***********************************************************************
*  expression tree and canonicalization tests
*************************************************************************/

void test_free_expr_no_crash() {
    free_expr(NULL);
    free_expr(create_const(42));
    free_expr(create_var("x"));
    free_expr(create_add(create_const(1), create_var("y")));
    free_expr(create_neg(create_var("x")));
    free_expr(create_mul(
        create_add(create_var("x"), create_const(2)),
        create_pow(create_var("x"), create_const(3))));
    printf("  test_free_expr_no_crash PASSED\n");
}

void test_canonicalize_constants() {
    // 3 + 4 => 7
    Expr* e = create_add(create_const(3), create_const(4));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);
    assert(c->poly->head && rat_eq(c->poly->head->coeff, R(7)));
    assert(c->poly->head->exponent == 0);
    free_expr(e);
    free_expr(c);
    printf("  test_canonicalize_constants PASSED\n");
}

void test_canonicalize_polynomial() {
    // x * x + 2 * x + 1 => x^2 + 2x + 1
    Expr* e = create_add(
        create_add(
            create_mul(create_var("x"), create_var("x")),
            create_mul(create_const(2), create_var("x"))),
        create_const(1));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);

    Term* t = c->poly->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(2)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    free_expr(e);
    free_expr(c);
    printf("  test_canonicalize_polynomial PASSED\n");
}

void test_canonicalize_negation() {
    // -(x + 2) => -x - 2
    Expr* e = create_neg(create_add(create_var("x"), create_const(2)));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);

    Term* t = c->poly->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(-1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(-2)));
    assert(t->next == NULL);

    free_expr(e);
    free_expr(c);
    printf("  test_canonicalize_negation PASSED\n");
}

void test_canonicalize_power() {
    // (x + 1)^2 => x^2 + 2x + 1
    Expr* e = create_pow(
        create_add(create_var("x"), create_const(1)),
        create_const(2));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);

    Term* t = c->poly->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(2)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));

    free_expr(e);
    free_expr(c);
    printf("  test_canonicalize_power PASSED\n");
}

void test_canonicalize_division() {
    // (x^2 - 1) / (x - 1) should reduce to x + 1
    Expr* num = create_sub(
        create_mul(create_var("x"), create_var("x")),
        create_const(1));
    Expr* den = create_sub(create_var("x"), create_const(1));
    Expr* e = create_div(num, den);
    Expr* c = canonicalize(e);

    // should reduce to polynomial x + 1
    assert(c && c->type == EXPR_POLY);
    Term* t = c->poly->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));

    free_expr(e);
    free_expr(c);
    printf("  test_canonicalize_division PASSED\n");
}

void test_simplify_combines_like_terms() {
    // 5x^3 + 3x^3 => 8x^3
    Expr* e = create_add(
        create_mul(create_const(5), create_pow(create_var("x"), create_const(3))),
        create_mul(create_const(3), create_pow(create_var("x"), create_const(3))));
    Expr* s = simplify(e);
    assert(s && s->type == EXPR_POLY);

    Term* t = s->poly->head;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(8)));
    assert(t->next == NULL);

    free_expr(e);
    free_expr(s);
    printf("  test_simplify_combines_like_terms PASSED\n");
}

/***********************************************************************
*  run all tests
*************************************************************************/

/***********************************************************************
*  factoring tests
*************************************************************************/

void test_factor_linear() {
    // x^2 - 5x + 6 = (x - 2)(x - 3)
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 2);
    add_term(p, R(-5), 1);
    add_term(p, R(6), 0);

    Factorization* f = factorize(p);
    assert(f && f->count == 2);
    assert(rat_is_one(f->content));

    // verify roots: factors should be (x-2) and (x-3) in some order
    // each factor is degree 1, multiplicity 1
    for (int i = 0; i < 2; i++) {
        assert(f->multiplicities[i] == 1);
        assert(polynomial_degree(f->factors[i]) == 1);
    }

    // verify by multiplying factors back: should give x^2 - 5x + 6
    Polynomial* product = poly_multiply(f->factors[0], f->factors[1]);
    Term* t = product->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(-5)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(6)));

    free_polynomial(product);
    free_factorization(f);
    free_polynomial(p);
    printf("  test_factor_linear PASSED\n");
}

void test_factor_repeated() {
    // x^3 - 3x^2 + 3x - 1 = (x - 1)^3
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 3);
    add_term(p, R(-3), 2);
    add_term(p, R(3), 1);
    add_term(p, R(-1), 0);

    Factorization* f = factorize(p);
    assert(f && f->count == 1);
    assert(f->multiplicities[0] == 3);

    // factor should be (x - 1)
    Term* t = f->factors[0]->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(-1)));

    free_factorization(f);
    free_polynomial(p);
    printf("  test_factor_repeated PASSED\n");
}

void test_factor_with_content() {
    // 6x^2 + 12x + 6 = 6(x + 1)^2
    Polynomial* p = create_polynomial();
    add_term(p, R(6), 2);
    add_term(p, R(12), 1);
    add_term(p, R(6), 0);

    Factorization* f = factorize(p);
    assert(f);
    assert(rat_eq(f->content, R(6)));
    assert(f->count == 1);
    assert(f->multiplicities[0] == 2);

    // factor should be (x + 1)
    Term* t = f->factors[0]->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));

    free_factorization(f);
    free_polynomial(p);
    printf("  test_factor_with_content PASSED\n");
}

void test_factor_irreducible() {
    // x^2 + 1 has no rational roots, should remain as single irreducible factor
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 2);
    add_term(p, R(1), 0);

    Factorization* f = factorize(p);
    assert(f && f->count == 1);
    assert(f->multiplicities[0] == 1);
    assert(polynomial_degree(f->factors[0]) == 2);

    free_factorization(f);
    free_polynomial(p);
    printf("  test_factor_irreducible PASSED\n");
}

void test_factor_rational_roots() {
    // 2x^3 - 3x^2 - 8x - 3 has roots 3, -1/2, -1
    // = 2(x - 3)(x + 1/2)(x + 1) = 2(x-3)(2x+1)(x+1)/2 ...
    // actually: 2(x-3)(x+1)(x+1/2) but let's verify via factorize
    Polynomial* p = create_polynomial();
    add_term(p, R(2), 3);
    add_term(p, R(-3), 2);
    add_term(p, R(-8), 1);
    add_term(p, R(-3), 0);

    Factorization* f = factorize(p);
    assert(f);
    // should have 3 linear factors
    int total_degree = 0;
    for (int i = 0; i < f->count; i++)
        total_degree += polynomial_degree(f->factors[i]) * f->multiplicities[i];
    assert(total_degree == 3);

    free_factorization(f);
    free_polynomial(p);
    printf("  test_factor_rational_roots PASSED\n");
}

void test_factor_to_string() {
    // x^2 - 5x + 6 = (x - 2)(x - 3)
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 2);
    add_term(p, R(-5), 1);
    add_term(p, R(6), 0);

    Factorization* f = factorize(p);
    char* s = factorization_to_string(f);
    // should contain parenthesized factors
    assert(strstr(s, "(") != NULL);
    printf("  factored form: %s\n", s);

    free(s);
    free_factorization(f);
    free_polynomial(p);
    printf("  test_factor_to_string PASSED\n");
}

/***********************************************************************
*  solver tests
*************************************************************************/

void test_solve_linear() {
    // 2x + 4 = 0 => x = -2
    Polynomial* p = create_polynomial();
    add_term(p, R(2), 1);
    add_term(p, R(4), 0);

    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 1);
    assert(rat_eq(sol->rational_roots[0], R(-2)));

    free_solutions(sol);
    free_polynomial(p);
    printf("  test_solve_linear PASSED\n");
}

void test_solve_quadratic_rational() {
    // x^2 - 5x + 6 = 0 => x = 2, x = 3
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 2);
    add_term(p, R(-5), 1);
    add_term(p, R(6), 0);

    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 2);

    // roots should be 2 and 3 (in some order)
    bool has_2 = false, has_3 = false;
    for (int i = 0; i < sol->num_rational; i++) {
        if (rat_eq(sol->rational_roots[i], R(2))) has_2 = true;
        if (rat_eq(sol->rational_roots[i], R(3))) has_3 = true;
    }
    assert(has_2 && has_3);

    free_solutions(sol);
    free_polynomial(p);
    printf("  test_solve_quadratic_rational PASSED\n");
}

void test_solve_quadratic_irrational() {
    // x^2 - 2 = 0 => x ~ ±1.414
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 2);
    add_term(p, R(-2), 0);

    Solutions* sol = solve_polynomial(p);
    // no rational roots, should have 2 irrational
    assert(sol->num_rational == 0);
    assert(sol->num_irrational == 2);
    // verify values
    double r1 = sol->irrational_roots[0];
    double r2 = sol->irrational_roots[1];
    assert(fabs(fabs(r1) - sqrt(2.0)) < 1e-10);
    assert(fabs(fabs(r2) - sqrt(2.0)) < 1e-10);
    assert(r1 * r2 < 0); // opposite signs

    free_solutions(sol);
    free_polynomial(p);
    printf("  test_solve_quadratic_irrational PASSED\n");
}

void test_solve_cubic() {
    // x^3 - 6x^2 + 11x - 6 = 0 => x = 1, 2, 3
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 3);
    add_term(p, R(-6), 2);
    add_term(p, R(11), 1);
    add_term(p, R(-6), 0);

    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 3);

    bool has_1 = false, has_2 = false, has_3 = false;
    for (int i = 0; i < sol->num_rational; i++) {
        if (rat_eq(sol->rational_roots[i], R(1))) has_1 = true;
        if (rat_eq(sol->rational_roots[i], R(2))) has_2 = true;
        if (rat_eq(sol->rational_roots[i], R(3))) has_3 = true;
    }
    assert(has_1 && has_2 && has_3);

    char* s = solutions_to_string(sol);
    printf("  solved: %s\n", s);
    free(s);

    free_solutions(sol);
    free_polynomial(p);
    printf("  test_solve_cubic PASSED\n");
}

void test_solve_no_real_roots() {
    // x^2 + 1 = 0 => no real roots
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 2);
    add_term(p, R(1), 0);

    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 0);
    assert(sol->num_irrational == 0);
    assert(sol->num_irreducible == 1);

    free_solutions(sol);
    free_polynomial(p);
    printf("  test_solve_no_real_roots PASSED\n");
}

/***********************************************************************
*  calculus tests
*************************************************************************/

void test_differentiate_poly() {
    // d/dx(3x^2 + 4x + 2) = 6x + 4 via expression-level differentiate
    Expr* e = create_add(
        create_add(
            create_mul(create_const(3), create_pow(create_var("x"), create_const(2))),
            create_mul(create_const(4), create_var("x"))),
        create_const(2));
    Expr* d = differentiate(e, "x");
    assert(d && d->type == EXPR_POLY);

    Term* t = d->poly->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(6)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(4)));
    assert(t->next == NULL);

    free_expr(e);
    free_expr(d);
    printf("  test_differentiate_poly PASSED\n");
}

void test_differentiate_via_parser() {
    // parse "x^3 + 2x^2 + x", differentiate => 3x^2 + 4x + 1
    Expr* parsed = parse_and_simplify("x^3 + 2x^2 + x");
    assert(parsed);
    Expr* d = differentiate(parsed, "x");
    assert(d && d->type == EXPR_POLY);

    Term* t = d->poly->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(3)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(4)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    free_expr(parsed);
    free_expr(d);
    printf("  test_differentiate_via_parser PASSED\n");
}

void test_integrate_poly() {
    // integrate 3x^2 + 2x => x^3 + x^2
    Expr* e = create_add(
        create_mul(create_const(3), create_pow(create_var("x"), create_const(2))),
        create_mul(create_const(2), create_var("x")));
    Expr* i = integrate(e, "x");
    assert(i && i->type == EXPR_POLY);

    Term* t = i->poly->head;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    free_expr(e);
    free_expr(i);
    printf("  test_integrate_poly PASSED\n");
}

void run_tests() {
    printf("Rational Number Tests:\n");
    test_rat_normalization();
    test_rat_arithmetic();
    test_rat_pow();
    test_rat_comparison();
    test_rat_to_string();
    test_rat_to_double();

    printf("\nPolynomial Arithmetic Tests:\n");
    test_polynomial_addition();
    test_polynomial_subtraction();
    test_polynomial_multiplication();
    test_polynomial_division();
    test_polynomial_divmod();
    test_polynomial_gcd();
    test_polynomial_derivative();
    test_polynomial_integral();
    test_polynomial_evaluate();
    test_polynomial_degree();
    test_poly_to_string();
    test_poly_copy();

    printf("\nParser Tests:\n");
    test_parse_implicit_mul();
    test_parse_unary_minus();
    test_parse_and_simplify_polynomial();

    printf("\nExpression Tree & Canonicalization Tests:\n");
    test_free_expr_no_crash();
    test_canonicalize_constants();
    test_canonicalize_polynomial();
    test_canonicalize_negation();
    test_canonicalize_power();
    test_canonicalize_division();
    test_simplify_combines_like_terms();

    printf("\nFactoring Tests:\n");
    test_factor_linear();
    test_factor_repeated();
    test_factor_with_content();
    test_factor_irreducible();
    test_factor_rational_roots();
    test_factor_to_string();

    printf("\nSolver Tests:\n");
    test_solve_linear();
    test_solve_quadratic_rational();
    test_solve_quadratic_irrational();
    test_solve_cubic();
    test_solve_no_real_roots();

    printf("\nCalculus Tests:\n");
    test_differentiate_poly();
    test_differentiate_via_parser();
    test_integrate_poly();

    printf("\nAll tests passed!\n");
}

int main() {
    run_tests();
    return 0;
}
