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
#include "matrix.h"
#include "partial_frac.h"
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

/***********************************************************************
*  complex / stress tests
*************************************************************************/

void test_parse_expand_product() {
    // (x + 1)(x + 2)(x + 3) via parser => x^3 + 6x^2 + 11x + 6
    Expr* e = parse_and_simplify("(x + 1) * (x + 2) * (x + 3)");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(6)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(11)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(6)));
    assert(t->next == NULL);
    free_expr(e);
    printf("  test_parse_expand_product PASSED\n");
}

void test_expand_then_factor_roundtrip() {
    // expand (x-1)(x-2)(x-3)(x-4) then factor it back
    // x^4 - 10x^3 + 35x^2 - 50x + 24
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 4);
    add_term(p, R(-10), 3);
    add_term(p, R(35), 2);
    add_term(p, R(-50), 1);
    add_term(p, R(24), 0);

    Factorization* f = factorize(p);
    assert(f && f->count == 4);
    for (int i = 0; i < 4; i++) {
        assert(f->multiplicities[i] == 1);
        assert(polynomial_degree(f->factors[i]) == 1);
    }

    // multiply all factors back together and verify
    Polynomial* rebuilt = create_polynomial();
    add_term(rebuilt, rat_one(), 0);
    for (int i = 0; i < f->count; i++) {
        Polynomial* tmp = poly_multiply(rebuilt, f->factors[i]);
        free_polynomial(rebuilt);
        rebuilt = tmp;
    }
    // should match original
    Term* t = rebuilt->head;
    assert(t && t->exponent == 4 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(-10)));
    t = t->next;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(35)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(-50)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(24)));

    free_polynomial(rebuilt);
    free_factorization(f);
    free_polynomial(p);
    printf("  test_expand_then_factor_roundtrip PASSED\n");
}

void test_differentiate_integrate_inverse() {
    // d/dx(integrate(f)) should give back f
    // f = 4x^3 - 6x^2 + 2x - 7
    Polynomial* f = create_polynomial();
    add_term(f, R(4), 3);
    add_term(f, R(-6), 2);
    add_term(f, R(2), 1);
    add_term(f, R(-7), 0);

    Polynomial* integral = poly_integral(f);
    Polynomial* deriv = poly_derivative(integral);

    // deriv should equal f exactly
    Term* tf = f->head;
    Term* td = deriv->head;
    while (tf && td) {
        assert(tf->exponent == td->exponent);
        assert(rat_eq(tf->coeff, td->coeff));
        tf = tf->next;
        td = td->next;
    }
    assert(tf == NULL && td == NULL);

    free_polynomial(f);
    free_polynomial(integral);
    free_polynomial(deriv);
    printf("  test_differentiate_integrate_inverse PASSED\n");
}

void test_high_degree_polynomial_multiply() {
    // (x^5 + 1) * (x^5 - 1) = x^10 - 1  (difference of squares)
    Polynomial* a = create_polynomial();
    add_term(a, R(1), 5);
    add_term(a, R(1), 0);

    Polynomial* b = create_polynomial();
    add_term(b, R(1), 5);
    add_term(b, R(-1), 0);

    Polynomial* product = poly_multiply(a, b);
    Term* t = product->head;
    assert(t && t->exponent == 10 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(-1)));
    assert(t->next == NULL);

    free_polynomial(a);
    free_polynomial(b);
    free_polynomial(product);
    printf("  test_high_degree_polynomial_multiply PASSED\n");
}

void test_rational_coefficient_factoring() {
    // (1/2)x^2 - (1/2) = (1/2)(x-1)(x+1)
    Polynomial* p = create_polynomial();
    add_term(p, rat_create(1, 2), 2);
    add_term(p, rat_create(-1, 2), 0);

    Factorization* f = factorize(p);
    assert(f);
    assert(rat_eq(f->content, rat_create(1, 2)));
    assert(f->count == 2);
    for (int i = 0; i < 2; i++) {
        assert(polynomial_degree(f->factors[i]) == 1);
        assert(f->multiplicities[i] == 1);
    }

    free_factorization(f);
    free_polynomial(p);
    printf("  test_rational_coefficient_factoring PASSED\n");
}

void test_solve_quartic() {
    // (x-1)(x+1)(x-2)(x+2) = x^4 - 5x^2 + 4
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 4);
    add_term(p, R(-5), 2);
    add_term(p, R(4), 0);

    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 4);

    bool has[5] = {false};
    for (int i = 0; i < sol->num_rational; i++) {
        rat_int_t v = sol->rational_roots[i].num;
        if (v >= -2 && v <= 2 && sol->rational_roots[i].den == 1)
            has[v + 2] = true;
    }
    assert(has[0] && has[1] && has[3] && has[4]); // -2,-1,1,2

    free_solutions(sol);
    free_polynomial(p);
    printf("  test_solve_quartic PASSED\n");
}

void test_solve_with_rational_root() {
    // 6x^3 - 11x^2 + 6x - 1 = 0  has roots 1/2, 1/3, 1
    // verify: (x-1)(2x-1)(3x-1) = 6x^3 - 11x^2 + 6x - 1
    Polynomial* p = create_polynomial();
    add_term(p, R(6), 3);
    add_term(p, R(-11), 2);
    add_term(p, R(6), 1);
    add_term(p, R(-1), 0);

    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 3);

    bool has_half = false, has_third = false, has_one = false;
    for (int i = 0; i < sol->num_rational; i++) {
        if (rat_eq(sol->rational_roots[i], rat_create(1, 2))) has_half = true;
        if (rat_eq(sol->rational_roots[i], rat_create(1, 3))) has_third = true;
        if (rat_eq(sol->rational_roots[i], R(1))) has_one = true;
    }
    assert(has_half && has_third && has_one);

    char* s = solutions_to_string(sol);
    printf("  solved: %s\n", s);
    free(s);

    free_solutions(sol);
    free_polynomial(p);
    printf("  test_solve_with_rational_root PASSED\n");
}

void test_gcd_nontrivial() {
    // gcd(x^4 - 1, x^6 - 1) = x^2 - 1
    // x^4-1 = (x^2-1)(x^2+1), x^6-1 = (x^2-1)(x^4+x^2+1)
    Polynomial* a = create_polynomial();
    add_term(a, R(1), 4);
    add_term(a, R(-1), 0);

    Polynomial* b = create_polynomial();
    add_term(b, R(1), 6);
    add_term(b, R(-1), 0);

    Polynomial* g = poly_gcd(a, b);

    // should be x^2 - 1 (monic)
    Term* t = g->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(-1)));
    assert(t->next == NULL);

    free_polynomial(a);
    free_polynomial(b);
    free_polynomial(g);
    printf("  test_gcd_nontrivial PASSED\n");
}

void test_definite_integral() {
    // integral of 3x^2 from 0 to 2 = [x^3]_0^2 = 8
    Polynomial* p = create_polynomial();
    add_term(p, R(3), 2);
    double result = poly_definite_integral(p, 0.0, 2.0);
    assert(fabs(result - 8.0) < 1e-10);
    free_polynomial(p);

    // integral of x^2 + x from 1 to 3 = [x^3/3 + x^2/2]_1^3
    // = (9 + 4.5) - (1/3 + 0.5) = 13.5 - 0.8333... = 12.6666...
    Polynomial* q = create_polynomial();
    add_term(q, R(1), 2);
    add_term(q, R(1), 1);
    double r2 = poly_definite_integral(q, 1.0, 3.0);
    assert(fabs(r2 - 38.0/3.0) < 1e-10);
    free_polynomial(q);

    printf("  test_definite_integral PASSED\n");
}

void test_sparse_polynomial_eval() {
    // x^100 + 1 evaluated at x=1 should be 2
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 100);
    add_term(p, R(1), 0);
    double val = poly_evaluate(p, 1.0);
    assert(fabs(val - 2.0) < 1e-10);

    // x^100 evaluated at x=0 should be 0
    double val0 = poly_evaluate(p, 0.0);
    assert(fabs(val0 - 1.0) < 1e-10); // constant term is 1

    free_polynomial(p);
    printf("  test_sparse_polynomial_eval PASSED\n");
}

void test_polynomial_chain_operations() {
    // compute (x+1)^4 by repeated multiplication, then differentiate twice
    // (x+1)^4 = x^4 + 4x^3 + 6x^2 + 4x + 1
    // d/dx = 4x^3 + 12x^2 + 12x + 4
    // d2/dx2 = 12x^2 + 24x + 12
    Polynomial* base = create_polynomial();
    add_term(base, R(1), 1);
    add_term(base, R(1), 0);

    Polynomial* p = poly_copy(base);
    for (int i = 0; i < 3; i++) {
        Polynomial* tmp = poly_multiply(p, base);
        free_polynomial(p);
        p = tmp;
    }

    // verify (x+1)^4
    Term* t = p->head;
    assert(t && t->exponent == 4 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(4)));
    t = t->next;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(6)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(4)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));

    Polynomial* d1 = poly_derivative(p);
    Polynomial* d2 = poly_derivative(d1);

    // 12x^2 + 24x + 12
    t = d2->head;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(12)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(24)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(12)));
    assert(t->next == NULL);

    free_polynomial(base);
    free_polynomial(p);
    free_polynomial(d1);
    free_polynomial(d2);
    printf("  test_polynomial_chain_operations PASSED\n");
}

void test_solve_degree5_partial() {
    // x^5 - x = x(x^4 - 1) = x(x^2-1)(x^2+1) = x(x-1)(x+1)(x^2+1)
    // rational roots: 0, 1, -1; irreducible: x^2+1
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 5);
    add_term(p, R(-1), 1);

    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 3);

    bool has_0 = false, has_1 = false, has_m1 = false;
    for (int i = 0; i < sol->num_rational; i++) {
        if (rat_eq(sol->rational_roots[i], R(0))) has_0 = true;
        if (rat_eq(sol->rational_roots[i], R(1))) has_1 = true;
        if (rat_eq(sol->rational_roots[i], R(-1))) has_m1 = true;
    }
    assert(has_0 && has_1 && has_m1);
    assert(sol->num_irreducible == 1);

    char* s = solutions_to_string(sol);
    printf("  solved: %s\n", s);
    free(s);

    free_solutions(sol);
    free_polynomial(p);
    printf("  test_solve_degree5_partial PASSED\n");
}

void test_parse_complex_expression() {
    // (2x^2 + 3x - 5) * (x - 1) - (x^3 + 2x^2 - 8x + 5)
    // = (2x^3 + 3x^2 - 5x - 2x^2 - 3x + 5) - (x^3 + 2x^2 - 8x + 5)
    // = (2x^3 + x^2 - 8x + 5) - (x^3 + 2x^2 - 8x + 5)
    // = x^3 - x^2
    Expr* e = parse_and_simplify("(2x^2 + 3x - 5) * (x - 1) - (x^3 + 2x^2 - 8x + 5)");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(-1)));
    assert(t->next == NULL);

    free_expr(e);
    printf("  test_parse_complex_expression PASSED\n");
}

void test_rational_arithmetic_stress() {
    // chain of operations that would accumulate error with doubles
    // ((1/7 + 1/11) * 77) should be exactly 18
    Rational a = rat_create(1, 7);
    Rational b = rat_create(1, 11);
    Rational sum = rat_add(a, b);            // 18/77
    Rational result = rat_mul(sum, rat_from_int(77)); // 18
    assert(rat_eq(result, R(18)));

    // (1/3)^10 * 3^10 = 1
    Rational third = rat_create(1, 3);
    Rational t10 = rat_pow_int(third, 10);   // 1/59049
    Rational three10 = rat_pow_int(R(3), 10); // 59049
    Rational product = rat_mul(t10, three10);
    assert(rat_eq(product, rat_one()));

    printf("  test_rational_arithmetic_stress PASSED\n");
}

void test_factor_degree6() {
    // x^6 - 1 = (x-1)(x+1)(x^4+x^2+1)
    // the quartic has no rational roots so stays as one irreducible factor
    Polynomial* p = create_polynomial();
    add_term(p, R(1), 6);
    add_term(p, R(-1), 0);

    Factorization* f = factorize(p);
    assert(f);

    int linear_count = 0;
    int total_degree = 0;
    for (int i = 0; i < f->count; i++) {
        int d = polynomial_degree(f->factors[i]);
        total_degree += d * f->multiplicities[i];
        if (d == 1) linear_count++;
    }
    assert(total_degree == 6);
    assert(linear_count == 2); // (x-1) and (x+1)

    // verify by evaluation: each root should zero the original
    assert(fabs(poly_evaluate(p, 1.0)) < 1e-10);
    assert(fabs(poly_evaluate(p, -1.0)) < 1e-10);

    // verify roundtrip: multiply all factors, scale by content
    Polynomial* rebuilt = create_polynomial();
    add_term(rebuilt, f->content, 0);
    for (int i = 0; i < f->count; i++) {
        for (int m = 0; m < f->multiplicities[i]; m++) {
            Polynomial* tmp = poly_multiply(rebuilt, f->factors[i]);
            free_polynomial(rebuilt);
            rebuilt = tmp;
        }
    }
    Term* t = rebuilt->head;
    assert(t && t->exponent == 6 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(-1)));
    assert(t->next == NULL);

    char* s = factorization_to_string(f);
    printf("  x^6-1 factored: %s\n", s);
    free(s);

    free_polynomial(rebuilt);
    free_factorization(f);
    free_polynomial(p);
    printf("  test_factor_degree6 PASSED\n");
}

void test_divmod_with_remainder() {
    // (x^3 + 1) / (x^2 + 1) = quotient x, remainder -x + 1
    Polynomial* a = create_polynomial();
    add_term(a, R(1), 3);
    add_term(a, R(1), 0);

    Polynomial* b = create_polynomial();
    add_term(b, R(1), 2);
    add_term(b, R(1), 0);

    PolyDivResult dr = poly_divmod(a, b);

    // quotient = x
    Term* t = dr.quotient->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    // remainder = -x + 1
    t = dr.remainder->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(-1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    // verify: a = b * quotient + remainder
    Polynomial* bq = poly_multiply(b, dr.quotient);
    Polynomial* check = poly_add(bq, dr.remainder);
    t = check->head;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(1)));
    assert(t->next == NULL);

    free_polynomial(a);
    free_polynomial(b);
    free_polynomial(bq);
    free_polynomial(check);
    free_polynomial(dr.quotient);
    free_polynomial(dr.remainder);
    printf("  test_divmod_with_remainder PASSED\n");
}

void test_canonicalize_nested_division() {
    // ((x^2 - 4) / (x - 2)) should simplify to x + 2
    // built from tree: (x*x - 4) / (x - 2)
    Expr* e = create_div(
        create_sub(create_pow(create_var("x"), create_const(2)), create_const(4)),
        create_sub(create_var("x"), create_const(2)));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);
    Term* t = c->poly->head;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(2)));
    assert(t->next == NULL);

    free_expr(e);
    free_expr(c);
    printf("  test_canonicalize_nested_division PASSED\n");
}

void test_power_of_binomial() {
    // (x + 2)^5 via canonicalization
    // = x^5 + 10x^4 + 40x^3 + 80x^2 + 80x + 32
    Expr* e = create_pow(
        create_add(create_var("x"), create_const(2)),
        create_const(5));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);

    Term* t = c->poly->head;
    assert(t && t->exponent == 5 && rat_eq(t->coeff, R(1)));
    t = t->next;
    assert(t && t->exponent == 4 && rat_eq(t->coeff, R(10)));
    t = t->next;
    assert(t && t->exponent == 3 && rat_eq(t->coeff, R(40)));
    t = t->next;
    assert(t && t->exponent == 2 && rat_eq(t->coeff, R(80)));
    t = t->next;
    assert(t && t->exponent == 1 && rat_eq(t->coeff, R(80)));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(t->coeff, R(32)));
    assert(t->next == NULL);

    free_expr(e);
    free_expr(c);
    printf("  test_power_of_binomial PASSED\n");
}

/***********************************************************************
*  matrix / linear systems tests
*************************************************************************/

void test_matrix_create_free() {
    Matrix* m = matrix_create(2, 3);
    assert(m->rows == 2 && m->cols == 3);
    // all entries should be zero
    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 3; c++)
            assert(rat_is_zero(MAT_AT(m, r, c)));
    // set and read back
    MAT_AT(m, 0, 1) = rat_from_int(7);
    assert(rat_eq(MAT_AT(m, 0, 1), R(7)));
    matrix_free(m);
    printf("  test_matrix_create_free PASSED\n");
}

void test_matrix_copy() {
    Matrix* m = matrix_create(2, 2);
    MAT_AT(m, 0, 0) = R(1); MAT_AT(m, 0, 1) = R(2);
    MAT_AT(m, 1, 0) = R(3); MAT_AT(m, 1, 1) = R(4);
    Matrix* c = matrix_copy(m);
    MAT_AT(m, 0, 0) = R(99); // modify original
    assert(rat_eq(MAT_AT(c, 0, 0), R(1))); // copy unaffected
    matrix_free(m);
    matrix_free(c);
    printf("  test_matrix_copy PASSED\n");
}

void test_vector_create_free() {
    Vector* v = vector_create(3);
    v->data[0] = R(5);
    v->data[1] = R(-3);
    v->data[2] = rat_create(1, 2);
    assert(rat_eq(v->data[2], rat_create(1, 2)));
    vector_free(v);
    printf("  test_vector_create_free PASSED\n");
}

void test_matrix_to_string() {
    Matrix* m = matrix_create(2, 2);
    MAT_AT(m, 0, 0) = R(1); MAT_AT(m, 0, 1) = R(2);
    MAT_AT(m, 1, 0) = R(3); MAT_AT(m, 1, 1) = R(4);
    char* s = matrix_to_string(m);
    assert(strstr(s, "1") != NULL);
    assert(strstr(s, "4") != NULL);
    free(s);
    matrix_free(m);
    printf("  test_matrix_to_string PASSED\n");
}

void test_rref_identity() {
    Matrix* m = matrix_create(3, 3);
    MAT_AT(m, 0, 0) = R(1);
    MAT_AT(m, 1, 1) = R(1);
    MAT_AT(m, 2, 2) = R(1);
    Matrix* rref = matrix_rref(m);
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            assert(rat_eq(MAT_AT(rref, r, c), (r == c) ? rat_one() : rat_zero()));
    matrix_free(m);
    matrix_free(rref);
    printf("  test_rref_identity PASSED\n");
}

void test_rref_2x3() {
    // [[1, 2, 3], [4, 5, 6]] → RREF [[1, 0, -1], [0, 1, 2]]
    Matrix* m = matrix_create(2, 3);
    MAT_AT(m, 0, 0) = R(1); MAT_AT(m, 0, 1) = R(2); MAT_AT(m, 0, 2) = R(3);
    MAT_AT(m, 1, 0) = R(4); MAT_AT(m, 1, 1) = R(5); MAT_AT(m, 1, 2) = R(6);
    Matrix* rref = matrix_rref(m);
    assert(rat_eq(MAT_AT(rref, 0, 0), R(1)));
    assert(rat_eq(MAT_AT(rref, 0, 1), R(0)));
    assert(rat_eq(MAT_AT(rref, 0, 2), R(-1)));
    assert(rat_eq(MAT_AT(rref, 1, 0), R(0)));
    assert(rat_eq(MAT_AT(rref, 1, 1), R(1)));
    assert(rat_eq(MAT_AT(rref, 1, 2), R(2)));
    matrix_free(m);
    matrix_free(rref);
    printf("  test_rref_2x3 PASSED\n");
}

void test_rref_needs_swap() {
    // [[0, 1], [1, 0]] requires row swap
    Matrix* m = matrix_create(2, 2);
    MAT_AT(m, 0, 0) = R(0); MAT_AT(m, 0, 1) = R(1);
    MAT_AT(m, 1, 0) = R(1); MAT_AT(m, 1, 1) = R(0);
    Matrix* rref = matrix_rref(m);
    assert(rat_eq(MAT_AT(rref, 0, 0), R(1)));
    assert(rat_eq(MAT_AT(rref, 1, 1), R(1)));
    matrix_free(m);
    matrix_free(rref);
    printf("  test_rref_needs_swap PASSED\n");
}

void test_rank() {
    // rank of [[1,2,3],[4,5,6],[7,8,9]] = 2 (rows are linearly dependent)
    Matrix* m = matrix_create(3, 3);
    MAT_AT(m, 0, 0) = R(1); MAT_AT(m, 0, 1) = R(2); MAT_AT(m, 0, 2) = R(3);
    MAT_AT(m, 1, 0) = R(4); MAT_AT(m, 1, 1) = R(5); MAT_AT(m, 1, 2) = R(6);
    MAT_AT(m, 2, 0) = R(7); MAT_AT(m, 2, 1) = R(8); MAT_AT(m, 2, 2) = R(9);
    assert(matrix_rank(m) == 2);
    matrix_free(m);
    printf("  test_rank PASSED\n");
}

void test_det_2x2() {
    // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
    Matrix* m = matrix_create(2, 2);
    MAT_AT(m, 0, 0) = R(1); MAT_AT(m, 0, 1) = R(2);
    MAT_AT(m, 1, 0) = R(3); MAT_AT(m, 1, 1) = R(4);
    Rational d = matrix_det(m);
    assert(rat_eq(d, R(-2)));
    matrix_free(m);
    printf("  test_det_2x2 PASSED\n");
}

void test_det_3x3() {
    // det([[2,1,1],[1,3,2],[1,0,0]]) = 2(0-0) - 1(0-2) + 1(0-3) = 0+2-3 = -1
    Matrix* m = matrix_create(3, 3);
    MAT_AT(m, 0, 0) = R(2); MAT_AT(m, 0, 1) = R(1); MAT_AT(m, 0, 2) = R(1);
    MAT_AT(m, 1, 0) = R(1); MAT_AT(m, 1, 1) = R(3); MAT_AT(m, 1, 2) = R(2);
    MAT_AT(m, 2, 0) = R(1); MAT_AT(m, 2, 1) = R(0); MAT_AT(m, 2, 2) = R(0);
    Rational d = matrix_det(m);
    assert(rat_eq(d, R(-1)));
    matrix_free(m);
    printf("  test_det_3x3 PASSED\n");
}

void test_det_singular() {
    // [[1,2],[2,4]] is singular
    Matrix* m = matrix_create(2, 2);
    MAT_AT(m, 0, 0) = R(1); MAT_AT(m, 0, 1) = R(2);
    MAT_AT(m, 1, 0) = R(2); MAT_AT(m, 1, 1) = R(4);
    assert(rat_is_zero(matrix_det(m)));
    matrix_free(m);
    printf("  test_det_singular PASSED\n");
}

void test_det_with_swap() {
    // [[0,1],[1,0]] needs swap, det = -1
    Matrix* m = matrix_create(2, 2);
    MAT_AT(m, 0, 0) = R(0); MAT_AT(m, 0, 1) = R(1);
    MAT_AT(m, 1, 0) = R(1); MAT_AT(m, 1, 1) = R(0);
    Rational d = matrix_det(m);
    assert(rat_eq(d, R(-1)));
    matrix_free(m);
    printf("  test_det_with_swap PASSED\n");
}

void test_solve_2x2_unique() {
    // x + 2y = 5, 3x + 4y = 11 => x = 1, y = 2
    Matrix* A = matrix_create(2, 2);
    MAT_AT(A, 0, 0) = R(1); MAT_AT(A, 0, 1) = R(2);
    MAT_AT(A, 1, 0) = R(3); MAT_AT(A, 1, 1) = R(4);
    Vector* b = vector_create(2);
    b->data[0] = R(5);
    b->data[1] = R(11);

    LinSysResult* res = linsys_solve(A, b);
    assert(res->status == LINSYS_UNIQUE);
    assert(rat_eq(res->solution->data[0], R(1)));
    assert(rat_eq(res->solution->data[1], R(2)));

    linsys_result_free(res);
    matrix_free(A);
    vector_free(b);
    printf("  test_solve_2x2_unique PASSED\n");
}

void test_solve_3x3_unique() {
    // x + y + z = 6, 2y + 5z = -4, 2x + 5y - z = 27
    // solution: x = 5, y = 3, z = -2
    Matrix* A = matrix_create(3, 3);
    MAT_AT(A, 0, 0) = R(1); MAT_AT(A, 0, 1) = R(1); MAT_AT(A, 0, 2) = R(1);
    MAT_AT(A, 1, 0) = R(0); MAT_AT(A, 1, 1) = R(2); MAT_AT(A, 1, 2) = R(5);
    MAT_AT(A, 2, 0) = R(2); MAT_AT(A, 2, 1) = R(5); MAT_AT(A, 2, 2) = R(-1);
    Vector* b = vector_create(3);
    b->data[0] = R(6); b->data[1] = R(-4); b->data[2] = R(27);

    LinSysResult* res = linsys_solve(A, b);
    assert(res->status == LINSYS_UNIQUE);
    assert(rat_eq(res->solution->data[0], R(5)));
    assert(rat_eq(res->solution->data[1], R(3)));
    assert(rat_eq(res->solution->data[2], R(-2)));

    linsys_result_free(res);
    matrix_free(A);
    vector_free(b);
    printf("  test_solve_3x3_unique PASSED\n");
}

void test_solve_inconsistent() {
    // x + y = 1, x + y = 2 — no solution
    Matrix* A = matrix_create(2, 2);
    MAT_AT(A, 0, 0) = R(1); MAT_AT(A, 0, 1) = R(1);
    MAT_AT(A, 1, 0) = R(1); MAT_AT(A, 1, 1) = R(1);
    Vector* b = vector_create(2);
    b->data[0] = R(1); b->data[1] = R(2);

    LinSysResult* res = linsys_solve(A, b);
    assert(res->status == LINSYS_INCONSISTENT);
    assert(res->solution == NULL);

    linsys_result_free(res);
    matrix_free(A);
    vector_free(b);
    printf("  test_solve_inconsistent PASSED\n");
}

void test_solve_infinite() {
    // x + 2y + 3z = 1 (1 equation, 3 unknowns)
    Matrix* A = matrix_create(1, 3);
    MAT_AT(A, 0, 0) = R(1); MAT_AT(A, 0, 1) = R(2); MAT_AT(A, 0, 2) = R(3);
    Vector* b = vector_create(1);
    b->data[0] = R(1);

    LinSysResult* res = linsys_solve(A, b);
    assert(res->status == LINSYS_INFINITE);
    assert(res->rank == 1);

    linsys_result_free(res);
    matrix_free(A);
    vector_free(b);
    printf("  test_solve_infinite PASSED\n");
}

void test_solve_hilbert_3x3() {
    // Hilbert matrix: H[i][j] = 1/(i+j+1), b = [1,1,1]
    // exact solution: x = [3, -24, 30]
    // floating-point solvers get this wrong (condition number ~524)
    Matrix* H = matrix_create(3, 3);
    MAT_AT(H, 0, 0) = rat_one();         MAT_AT(H, 0, 1) = rat_create(1, 2); MAT_AT(H, 0, 2) = rat_create(1, 3);
    MAT_AT(H, 1, 0) = rat_create(1, 2);  MAT_AT(H, 1, 1) = rat_create(1, 3); MAT_AT(H, 1, 2) = rat_create(1, 4);
    MAT_AT(H, 2, 0) = rat_create(1, 3);  MAT_AT(H, 2, 1) = rat_create(1, 4); MAT_AT(H, 2, 2) = rat_create(1, 5);

    Vector* b = vector_create(3);
    b->data[0] = rat_one(); b->data[1] = rat_one(); b->data[2] = rat_one();

    LinSysResult* res = linsys_solve(H, b);
    assert(res->status == LINSYS_UNIQUE);
    assert(rat_eq(res->solution->data[0], R(3)));
    assert(rat_eq(res->solution->data[1], R(-24)));
    assert(rat_eq(res->solution->data[2], R(30)));

    char* s = linsys_result_to_string(res);
    printf("  hilbert solved: %s\n", s);
    free(s);

    linsys_result_free(res);
    matrix_free(H);
    vector_free(b);
    printf("  test_solve_hilbert_3x3 PASSED\n");
}

void test_inverse_2x2() {
    // inverse of [[1,2],[3,4]], verify A * A^-1 = I
    Matrix* A = matrix_create(2, 2);
    MAT_AT(A, 0, 0) = R(1); MAT_AT(A, 0, 1) = R(2);
    MAT_AT(A, 1, 0) = R(3); MAT_AT(A, 1, 1) = R(4);

    Matrix* inv = matrix_inverse(A);
    assert(inv != NULL);

    // expected: [[-2, 1], [3/2, -1/2]]
    assert(rat_eq(MAT_AT(inv, 0, 0), R(-2)));
    assert(rat_eq(MAT_AT(inv, 0, 1), R(1)));
    assert(rat_eq(MAT_AT(inv, 1, 0), rat_create(3, 2)));
    assert(rat_eq(MAT_AT(inv, 1, 1), rat_create(-1, 2)));

    matrix_free(A);
    matrix_free(inv);
    printf("  test_inverse_2x2 PASSED\n");
}

void test_inverse_singular() {
    Matrix* m = matrix_create(2, 2);
    MAT_AT(m, 0, 0) = R(1); MAT_AT(m, 0, 1) = R(2);
    MAT_AT(m, 1, 0) = R(2); MAT_AT(m, 1, 1) = R(4);
    assert(matrix_inverse(m) == NULL);
    matrix_free(m);
    printf("  test_inverse_singular PASSED\n");
}

void test_matrix_parse() {
    Matrix* m = matrix_parse("[[1, 2], [3, 4]]");
    assert(m && m->rows == 2 && m->cols == 2);
    assert(rat_eq(MAT_AT(m, 0, 0), R(1)));
    assert(rat_eq(MAT_AT(m, 1, 1), R(4)));
    matrix_free(m);

    // rational entries
    Matrix* m2 = matrix_parse("[[1/2, -3/4], [0, 5]]");
    assert(m2 && m2->rows == 2 && m2->cols == 2);
    assert(rat_eq(MAT_AT(m2, 0, 0), rat_create(1, 2)));
    assert(rat_eq(MAT_AT(m2, 0, 1), rat_create(-3, 4)));
    matrix_free(m2);

    printf("  test_matrix_parse PASSED\n");
}

void test_vector_parse() {
    Vector* v = vector_parse("[1, -2, 3/5]");
    assert(v && v->n == 3);
    assert(rat_eq(v->data[0], R(1)));
    assert(rat_eq(v->data[1], R(-2)));
    assert(rat_eq(v->data[2], rat_create(3, 5)));
    vector_free(v);

    assert(vector_parse("not a vector") == NULL);
    printf("  test_vector_parse PASSED\n");
}

// helper: get coefficient at given degree from a polynomial
static Rational test_get_coeff(const Polynomial* poly, int degree) {
    Term* t = poly->head;
    while (t) {
        if (t->exponent == degree) return t->coeff;
        t = t->next;
    }
    return rat_zero();
}

/***********************************************************************
*  partial fraction tests
*************************************************************************/

void test_pf_simple() {
    // (2x+3) / (x^2-1) = (2x+3)/((x-1)(x+1))
    // = 5/2 * 1/(x-1) + (-1/2) * 1/(x+1)
    Polynomial* num = create_polynomial();
    add_term(num, R(2), 1);
    add_term(num, R(3), 0);

    Polynomial* den = create_polynomial();
    add_term(den, R(1), 2);
    add_term(den, R(-1), 0);

    PartialFractionDecomp* pf = partial_fractions(num, den);
    assert(pf);
    assert(pf->poly_part == NULL || poly_is_zero(pf->poly_part));
    assert(pf->num_terms == 2);

    // verify by checking that the decomposition is correct:
    // reconstruct and compare. Each term has coeffs for powers 1..m
    // For simple roots, each term has 1 coefficient
    char* s = pf_to_string(pf);
    printf("  (2x+3)/(x^2-1) = %s\n", s);
    free(s);

    // verify coefficients: find term for (x-1) and (x+1)
    for (int i = 0; i < pf->num_terms; i++) {
        assert(!pf->terms[i].is_quadratic);
        assert(pf->terms[i].power == 1);
        // check coeff: for (x-1), should be 5/2; for (x+1), should be -1/2
        Rational root_neg = test_get_coeff(pf->terms[i].factor, 0);
        if (rat_eq(root_neg, R(-1))) {
            // factor is (x - 1)
            assert(rat_eq(pf->terms[i].coeffs[0], rat_create(5, 2)));
        } else {
            // factor is (x + 1)
            assert(rat_eq(pf->terms[i].coeffs[0], rat_create(-1, 2)));
        }
    }

    free_pf_decomp(pf);
    free_polynomial(num);
    free_polynomial(den);
    printf("  test_pf_simple PASSED\n");
}

void test_pf_repeated_root() {
    // (3x+5) / (x-1)^2 = A/(x-1) + B/(x-1)^2
    // clear: 3x+5 = A(x-1) + B = Ax + (B-A)
    // A = 3, B-A = 5, so B = 8
    Polynomial* num = create_polynomial();
    add_term(num, R(3), 1);
    add_term(num, R(5), 0);

    Polynomial* den = create_polynomial();
    add_term(den, R(1), 2);
    add_term(den, R(-2), 1);
    add_term(den, R(1), 0);

    PartialFractionDecomp* pf = partial_fractions(num, den);
    assert(pf);
    assert(pf->num_terms == 1);
    assert(pf->terms[0].power == 2);
    assert(rat_eq(pf->terms[0].coeffs[0], R(3)));  // A₁ = 3 (for 1/(x-1))
    assert(rat_eq(pf->terms[0].coeffs[1], R(8)));  // A₂ = 8 (for 1/(x-1)²)

    char* s = pf_to_string(pf);
    printf("  (3x+5)/(x-1)^2 = %s\n", s);
    free(s);

    free_pf_decomp(pf);
    free_polynomial(num);
    free_polynomial(den);
    printf("  test_pf_repeated_root PASSED\n");
}

void test_pf_improper() {
    // (x^2+1) / (x^2-1) = 1 + 2/(x^2-1) = 1 + 1/(x-1) - 1/(x+1)
    Polynomial* num = create_polynomial();
    add_term(num, R(1), 2);
    add_term(num, R(1), 0);

    Polynomial* den = create_polynomial();
    add_term(den, R(1), 2);
    add_term(den, R(-1), 0);

    PartialFractionDecomp* pf = partial_fractions(num, den);
    assert(pf);
    // should have polynomial part = 1
    assert(pf->poly_part != NULL);
    assert(pf->poly_part->head != NULL);
    assert(rat_eq(pf->poly_part->head->coeff, R(1)));
    assert(pf->poly_part->head->exponent == 0);
    // should have 2 PF terms
    assert(pf->num_terms == 2);

    char* s = pf_to_string(pf);
    printf("  (x^2+1)/(x^2-1) = %s\n", s);
    free(s);

    free_pf_decomp(pf);
    free_polynomial(num);
    free_polynomial(den);
    printf("  test_pf_improper PASSED\n");
}

/***********************************************************************
*  rational function integration tests
*************************************************************************/

void test_integrate_rational_simple() {
    // ∫ (2x+3)/(x^2-1) dx = 5/2 ln|x-1| - 1/2 ln|x+1|
    Polynomial* num = create_polynomial();
    add_term(num, R(2), 1);
    add_term(num, R(3), 0);
    Polynomial* den = create_polynomial();
    add_term(den, R(1), 2);
    add_term(den, R(-1), 0);

    Expr* ratfn = create_rational_fn(num, den);
    Expr* result = integrate(ratfn, "x");
    assert(result != NULL);

    char* s = expr_to_string(result);
    printf("  ∫(2x+3)/(x²-1)dx = %s\n", s);
    // should contain "ln" terms
    assert(strstr(s, "ln") != NULL);
    free(s);

    free_expr(ratfn);
    free_expr(result);
    printf("  test_integrate_rational_simple PASSED\n");
}

void test_integrate_rational_repeated() {
    // ∫ 1/(x-1)^2 dx = -1/(x-1)
    Polynomial* num = create_polynomial();
    add_term(num, R(1), 0);
    Polynomial* den = create_polynomial();
    add_term(den, R(1), 2);
    add_term(den, R(-2), 1);
    add_term(den, R(1), 0);

    Expr* ratfn = create_rational_fn(num, den);
    Expr* result = integrate(ratfn, "x");
    assert(result != NULL);

    char* s = expr_to_string(result);
    printf("  ∫1/(x-1)²dx = %s\n", s);
    // should NOT contain "ln" (power rule, not log)
    assert(strstr(s, "ln") == NULL);
    free(s);

    free_expr(ratfn);
    free_expr(result);
    printf("  test_integrate_rational_repeated PASSED\n");
}

void test_integrate_rational_improper() {
    // ∫ (x^2+1)/(x^2-1) dx = x + ln terms
    Polynomial* num = create_polynomial();
    add_term(num, R(1), 2);
    add_term(num, R(1), 0);
    Polynomial* den = create_polynomial();
    add_term(den, R(1), 2);
    add_term(den, R(-1), 0);

    Expr* ratfn = create_rational_fn(num, den);
    Expr* result = integrate(ratfn, "x");
    assert(result != NULL);

    char* s = expr_to_string(result);
    printf("  ∫(x²+1)/(x²-1)dx = %s\n", s);
    // should contain both polynomial and ln parts
    assert(strstr(s, "ln") != NULL);
    free(s);

    free_expr(ratfn);
    free_expr(result);
    printf("  test_integrate_rational_improper PASSED\n");
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

    printf("\nComplex / Stress Tests:\n");
    test_parse_expand_product();
    test_expand_then_factor_roundtrip();
    test_differentiate_integrate_inverse();
    test_high_degree_polynomial_multiply();
    test_rational_coefficient_factoring();
    test_solve_quartic();
    test_solve_with_rational_root();
    test_gcd_nontrivial();
    test_definite_integral();
    test_sparse_polynomial_eval();
    test_polynomial_chain_operations();
    test_solve_degree5_partial();
    test_parse_complex_expression();
    test_rational_arithmetic_stress();
    test_factor_degree6();
    test_divmod_with_remainder();
    test_canonicalize_nested_division();
    test_power_of_binomial();

    printf("\nMatrix / Linear Systems Tests:\n");
    test_matrix_create_free();
    test_matrix_copy();
    test_vector_create_free();
    test_matrix_to_string();
    test_rref_identity();
    test_rref_2x3();
    test_rref_needs_swap();
    test_rank();
    test_det_2x2();
    test_det_3x3();
    test_det_singular();
    test_det_with_swap();
    test_solve_2x2_unique();
    test_solve_3x3_unique();
    test_solve_inconsistent();
    test_solve_infinite();
    test_solve_hilbert_3x3();
    test_inverse_2x2();
    test_inverse_singular();
    test_matrix_parse();
    test_vector_parse();

    printf("\nPartial Fraction Tests:\n");
    test_pf_simple();
    test_pf_repeated_root();
    test_pf_improper();

    printf("\nRational Integration Tests:\n");
    test_integrate_rational_simple();
    test_integrate_rational_repeated();
    test_integrate_rational_improper();

    printf("\nAll tests passed!\n");
}

int main() {
    run_tests();
    return 0;
}
