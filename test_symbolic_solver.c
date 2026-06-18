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

// helper macro: creates a stack-allocated Rational from int (auto-initialized)
#define R(n) rat_from_int(n)
// helper: assert coeff equals integer, then clear the temp
#define ASSERT_COEFF(term, expected_exp, expected_int) do { \
    Rational _tmp = R(expected_int); \
    assert((term) && (term)->exponent == (expected_exp) && rat_eq(&(term)->coeff, &_tmp)); \
    rat_clear(&_tmp); \
} while(0)

// helper: get coefficient at given degree from a polynomial
static void test_get_coeff(Rational* result, const Polynomial* poly, int degree) {
    Term* t = poly->head;
    while (t) {
        if (t->exponent == degree) { rat_set(result, &t->coeff); return; }
        t = t->next;
    }
    mpq_set_si(result->val, 0, 1);
}

/***********************************************************************
*  rational number tests
*************************************************************************/

void test_rat_normalization() {
    Rational r1 = rat_from_si(6, 4);
    assert(rat_num_si(&r1) == 3 && rat_den_si(&r1) == 2);

    Rational r2 = rat_from_si(-6, -4);
    assert(rat_num_si(&r2) == 3 && rat_den_si(&r2) == 2);

    Rational r3 = rat_from_si(6, -4);
    assert(rat_num_si(&r3) == -3 && rat_den_si(&r3) == 2);

    Rational r4 = rat_from_si(0, 5);
    assert(rat_is_zero(&r4));

    rat_clear(&r1); rat_clear(&r2); rat_clear(&r3); rat_clear(&r4);
    printf("  test_rat_normalization PASSED\n");
}

void test_rat_arithmetic() {
    Rational a = rat_from_si(1, 2), b = rat_from_si(1, 3);
    Rational sum; rat_init(&sum);
    rat_add(&sum, &a, &b);
    Rational expected = rat_from_si(5, 6);
    assert(rat_eq(&sum, &expected));
    rat_clear(&a); rat_clear(&b); rat_clear(&sum); rat_clear(&expected);

    a = rat_from_si(3, 4); b = rat_from_si(1, 4);
    Rational diff; rat_init(&diff);
    rat_sub(&diff, &a, &b);
    expected = rat_from_si(1, 2);
    assert(rat_eq(&diff, &expected));
    rat_clear(&a); rat_clear(&b); rat_clear(&diff); rat_clear(&expected);

    a = rat_from_si(2, 3); b = rat_from_si(3, 4);
    Rational prod; rat_init(&prod);
    rat_mul(&prod, &a, &b);
    expected = rat_from_si(1, 2);
    assert(rat_eq(&prod, &expected));
    rat_clear(&a); rat_clear(&b); rat_clear(&prod); rat_clear(&expected);

    a = rat_from_si(2, 3); b = rat_from_si(4, 5);
    Rational quot; rat_init(&quot);
    rat_div(&quot, &a, &b);
    expected = rat_from_si(5, 6);
    assert(rat_eq(&quot, &expected));
    rat_clear(&a); rat_clear(&b); rat_clear(&quot); rat_clear(&expected);

    a = rat_from_si(3, 7);
    Rational neg; rat_init(&neg);
    rat_neg(&neg, &a);
    expected = rat_from_si(-3, 7);
    assert(rat_eq(&neg, &expected));
    rat_clear(&a); rat_clear(&neg); rat_clear(&expected);

    printf("  test_rat_arithmetic PASSED\n");
}

void test_rat_pow() {
    Rational base = rat_from_si(2, 3);
    Rational r; rat_init(&r);
    rat_pow_int(&r, &base, 3);
    Rational expected = rat_from_si(8, 27);
    assert(rat_eq(&r, &expected));
    rat_clear(&expected);

    rat_pow_int(&r, &base, 0);
    expected = rat_one_val();
    assert(rat_eq(&r, &expected));
    rat_clear(&expected);

    rat_pow_int(&r, &base, -1);
    expected = rat_from_si(3, 2);
    assert(rat_eq(&r, &expected));

    rat_clear(&base); rat_clear(&r); rat_clear(&expected);
    printf("  test_rat_pow PASSED\n");
}

void test_rat_comparison() {
    Rational a = rat_from_si(2, 4), b = rat_from_si(1, 2);
    assert(rat_eq(&a, &b));
    rat_clear(&a); rat_clear(&b);

    a = rat_from_si(1, 2); b = rat_from_si(1, 3);
    assert(!rat_eq(&a, &b));
    rat_clear(&a); rat_clear(&b);

    a = rat_from_si(0, 42);
    assert(rat_is_zero(&a));
    rat_clear(&a);

    a = rat_from_si(5, 5);
    assert(rat_is_one(&a));
    rat_clear(&a);

    a = rat_from_si(-1, 3);
    assert(rat_is_negative(&a));
    rat_clear(&a);

    a = rat_from_si(1, 2); b = rat_from_si(2, 3);
    assert(rat_cmp(&a, &b) < 0);
    rat_clear(&a); rat_clear(&b);

    printf("  test_rat_comparison PASSED\n");
}

void test_rat_to_string() {
    Rational r = rat_from_si(3, 1);
    char* s = rat_to_string(&r);
    assert(strcmp(s, "3") == 0);
    free(s); rat_clear(&r);

    r = rat_from_si(3, 5);
    s = rat_to_string(&r);
    assert(strcmp(s, "3/5") == 0);
    free(s); rat_clear(&r);

    r = rat_from_si(-1, 2);
    s = rat_to_string(&r);
    assert(strcmp(s, "-1/2") == 0);
    free(s); rat_clear(&r);

    printf("  test_rat_to_string PASSED\n");
}

void test_rat_to_double() {
    Rational r = rat_from_si(1, 2);
    assert(fabs(rat_to_double(&r) - 0.5) < 1e-15);
    rat_clear(&r);

    r = rat_from_si(1, 3);
    assert(fabs(rat_to_double(&r) - 0.333333333333333) < 1e-10);
    rat_clear(&r);

    printf("  test_rat_to_double PASSED\n");
}

/***********************************************************************
*  polynomial arithmetic tests
*************************************************************************/

void test_polynomial_addition() {
    Polynomial* p1 = create_polynomial();
    Rational c; c = R(3); add_term(p1, &c, 2); rat_clear(&c);
    c = R(4); add_term(p1, &c, 1); rat_clear(&c);
    c = R(2); add_term(p1, &c, 0); rat_clear(&c);

    Polynomial* p2 = create_polynomial();
    c = R(5); add_term(p2, &c, 2); rat_clear(&c);
    c = R(3); add_term(p2, &c, 1); rat_clear(&c);
    c = R(1); add_term(p2, &c, 0); rat_clear(&c);

    Polynomial* sum = poly_add(p1, p2);
    Term* t = sum->head;
    ASSERT_COEFF(t, 2, 8); t = t->next;
    ASSERT_COEFF(t, 1, 7); t = t->next;
    ASSERT_COEFF(t, 0, 3); assert(t->next == NULL);

    free_polynomial(p1); free_polynomial(p2); free_polynomial(sum);
    printf("  test_polynomial_addition PASSED\n");
}

void test_polynomial_subtraction() {
    Polynomial* p1 = create_polynomial();
    Rational c; c = R(3); add_term(p1, &c, 2); rat_clear(&c);
    c = R(4); add_term(p1, &c, 1); rat_clear(&c);
    c = R(2); add_term(p1, &c, 0); rat_clear(&c);

    Polynomial* p2 = create_polynomial();
    c = R(5); add_term(p2, &c, 2); rat_clear(&c);
    c = R(3); add_term(p2, &c, 1); rat_clear(&c);
    c = R(1); add_term(p2, &c, 0); rat_clear(&c);

    Polynomial* diff = poly_subtract(p1, p2);
    Term* t = diff->head;
    ASSERT_COEFF(t, 2, -2); t = t->next;
    ASSERT_COEFF(t, 1, 1); t = t->next;
    ASSERT_COEFF(t, 0, 1); assert(t->next == NULL);

    free_polynomial(p1); free_polynomial(p2); free_polynomial(diff);
    printf("  test_polynomial_subtraction PASSED\n");
}

void test_polynomial_multiplication() {
    Polynomial* p1 = create_polynomial();
    Rational c; c = R(3); add_term(p1, &c, 2); rat_clear(&c);
    c = R(4); add_term(p1, &c, 1); rat_clear(&c);
    c = R(2); add_term(p1, &c, 0); rat_clear(&c);

    Polynomial* p2 = create_polynomial();
    c = R(5); add_term(p2, &c, 2); rat_clear(&c);
    c = R(3); add_term(p2, &c, 1); rat_clear(&c);
    c = R(1); add_term(p2, &c, 0); rat_clear(&c);

    Polynomial* product = poly_multiply(p1, p2);
    Term* t = product->head;
    ASSERT_COEFF(t, 4, 15); t = t->next;
    ASSERT_COEFF(t, 3, 29); t = t->next;
    ASSERT_COEFF(t, 2, 25); t = t->next;
    ASSERT_COEFF(t, 1, 10); t = t->next;
    ASSERT_COEFF(t, 0, 2); assert(t->next == NULL);

    free_polynomial(p1); free_polynomial(p2); free_polynomial(product);
    printf("  test_polynomial_multiplication PASSED\n");
}

void test_polynomial_division() {
    Polynomial* dividend = create_polynomial();
    Rational c; c = R(1); add_term(dividend, &c, 2); rat_clear(&c);
    c = R(-1); add_term(dividend, &c, 0); rat_clear(&c);

    Polynomial* divisor = create_polynomial();
    c = R(1); add_term(divisor, &c, 1); rat_clear(&c);
    c = R(-1); add_term(divisor, &c, 0); rat_clear(&c);

    Polynomial* quotient = poly_divide(dividend, divisor);
    Term* t = quotient->head;
    ASSERT_COEFF(t, 1, 1); t = t->next;
    ASSERT_COEFF(t, 0, 1); assert(t->next == NULL);

    free_polynomial(dividend); free_polynomial(divisor); free_polynomial(quotient);
    printf("  test_polynomial_division PASSED\n");
}

void test_polynomial_divmod() {
    Polynomial* dividend = create_polynomial();
    Rational c; c = R(1); add_term(dividend, &c, 3); rat_clear(&c);
    c = R(2); add_term(dividend, &c, 2); rat_clear(&c);
    c = R(3); add_term(dividend, &c, 1); rat_clear(&c);
    c = R(4); add_term(dividend, &c, 0); rat_clear(&c);

    Polynomial* divisor = create_polynomial();
    c = R(1); add_term(divisor, &c, 1); rat_clear(&c);
    c = R(1); add_term(divisor, &c, 0); rat_clear(&c);

    PolyDivResult dr = poly_divmod(dividend, divisor);
    Term* t = dr.quotient->head;
    ASSERT_COEFF(t, 2, 1); t = t->next;
    ASSERT_COEFF(t, 1, 1); t = t->next;
    ASSERT_COEFF(t, 0, 2); assert(t->next == NULL);

    t = dr.remainder->head;
    ASSERT_COEFF(t, 0, 2); assert(t->next == NULL);

    free_polynomial(dividend); free_polynomial(divisor);
    free_polynomial(dr.quotient); free_polynomial(dr.remainder);
    printf("  test_polynomial_divmod PASSED\n");
}

void test_polynomial_gcd() {
    Polynomial* p1 = create_polynomial();
    Rational c; c = R(1); add_term(p1, &c, 2); rat_clear(&c);
    c = R(-1); add_term(p1, &c, 0); rat_clear(&c);

    Polynomial* p2 = create_polynomial();
    c = R(1); add_term(p2, &c, 2); rat_clear(&c);
    c = R(-2); add_term(p2, &c, 1); rat_clear(&c);
    c = R(1); add_term(p2, &c, 0); rat_clear(&c);

    Polynomial* gcd = poly_gcd(p1, p2);
    Term* t = gcd->head;
    ASSERT_COEFF(t, 1, 1); t = t->next;
    ASSERT_COEFF(t, 0, -1); assert(t->next == NULL);

    free_polynomial(p1); free_polynomial(p2); free_polynomial(gcd);
    printf("  test_polynomial_gcd PASSED\n");
}

void test_polynomial_derivative() {
    Polynomial* poly = create_polynomial();
    Rational c; c = R(3); add_term(poly, &c, 2); rat_clear(&c);
    c = R(4); add_term(poly, &c, 1); rat_clear(&c);
    c = R(2); add_term(poly, &c, 0); rat_clear(&c);

    Polynomial* deriv = poly_derivative(poly);
    Term* t = deriv->head;
    ASSERT_COEFF(t, 1, 6); t = t->next;
    ASSERT_COEFF(t, 0, 4); assert(t->next == NULL);

    free_polynomial(poly); free_polynomial(deriv);
    printf("  test_polynomial_derivative PASSED\n");
}

void test_polynomial_integral() {
    Polynomial* poly = create_polynomial();
    Rational c; c = R(2); add_term(poly, &c, 1); rat_clear(&c);
    c = R(3); add_term(poly, &c, 0); rat_clear(&c);

    Polynomial* integral = poly_integral(poly);
    Term* t = integral->head;
    ASSERT_COEFF(t, 2, 1); t = t->next;
    ASSERT_COEFF(t, 1, 3); assert(t->next == NULL);

    free_polynomial(poly); free_polynomial(integral);
    printf("  test_polynomial_integral PASSED\n");
}

void test_polynomial_evaluate() {
    Polynomial* poly = create_polynomial();
    Rational c; c = R(3); add_term(poly, &c, 2); rat_clear(&c);
    c = R(4); add_term(poly, &c, 1); rat_clear(&c);
    c = R(2); add_term(poly, &c, 0); rat_clear(&c);

    assert(fabs(poly_evaluate(poly, 2.0) - 22.0) < 1e-10);
    free_polynomial(poly);
    printf("  test_polynomial_evaluate PASSED\n");
}

void test_polynomial_degree() {
    Polynomial* poly = create_polynomial();
    Rational c; c = R(3); add_term(poly, &c, 5); rat_clear(&c);
    c = R(1); add_term(poly, &c, 2); rat_clear(&c);
    c = R(7); add_term(poly, &c, 0); rat_clear(&c);
    assert(polynomial_degree(poly) == 5);
    free_polynomial(poly);
    printf("  test_polynomial_degree PASSED\n");
}

void test_poly_to_string() {
    Polynomial* poly = create_polynomial();
    Rational c; c = R(3); add_term(poly, &c, 2); rat_clear(&c);
    c = R(-4); add_term(poly, &c, 1); rat_clear(&c);
    c = R(5); add_term(poly, &c, 0); rat_clear(&c);
    char* s = poly_to_string(poly);
    assert(strcmp(s, "3x^2 - 4x + 5") == 0);
    free(s);

    Polynomial* p2 = create_polynomial();
    c = rat_from_si(1, 2); add_term(p2, &c, 2); rat_clear(&c);
    c = rat_from_si(-3, 4); add_term(p2, &c, 0); rat_clear(&c);
    char* s2 = poly_to_string(p2);
    assert(strcmp(s2, "1/2x^2 - 3/4") == 0);
    free(s2);

    free_polynomial(poly); free_polynomial(p2);
    printf("  test_poly_to_string PASSED\n");
}

void test_poly_copy() {
    Polynomial* poly = create_polynomial();
    Rational c; c = R(3); add_term(poly, &c, 2); rat_clear(&c);
    c = R(1); add_term(poly, &c, 0); rat_clear(&c);
    Polynomial* copy = poly_copy(poly);
    c = R(5); add_term(poly, &c, 1); rat_clear(&c);

    Term* t = copy->head;
    ASSERT_COEFF(t, 2, 3); t = t->next;
    ASSERT_COEFF(t, 0, 1); assert(t->next == NULL);

    free_polynomial(poly); free_polynomial(copy);
    printf("  test_poly_copy PASSED\n");
}

/***********************************************************************
*  parser tests
*************************************************************************/

void test_parse_implicit_mul() {
    Expr* e = parse_and_simplify("5x^3");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    ASSERT_COEFF(t, 3, 5); assert(t->next == NULL);
    free_expr(e);
    printf("  test_parse_implicit_mul PASSED\n");
}

void test_parse_unary_minus() {
    Expr* e = parse_and_simplify("-x");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    ASSERT_COEFF(t, 1, -1); assert(t->next == NULL);
    free_expr(e);
    printf("  test_parse_unary_minus PASSED\n");
}

void test_parse_and_simplify_polynomial() {
    Expr* e = parse_and_simplify("5x^3 + 3x^3 - 2x + 4 + 7");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    ASSERT_COEFF(t, 3, 8); t = t->next;
    ASSERT_COEFF(t, 1, -2); t = t->next;
    ASSERT_COEFF(t, 0, 11); assert(t->next == NULL);
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
    Expr* e = create_add(create_const(3), create_const(4));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);
    ASSERT_COEFF(c->poly->head, 0, 7);
    free_expr(e); free_expr(c);
    printf("  test_canonicalize_constants PASSED\n");
}

void test_canonicalize_polynomial() {
    Expr* e = create_add(
        create_add(create_mul(create_var("x"), create_var("x")),
                   create_mul(create_const(2), create_var("x"))),
        create_const(1));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);
    Term* t = c->poly->head;
    ASSERT_COEFF(t, 2, 1); t = t->next;
    ASSERT_COEFF(t, 1, 2); t = t->next;
    ASSERT_COEFF(t, 0, 1); assert(t->next == NULL);
    free_expr(e); free_expr(c);
    printf("  test_canonicalize_polynomial PASSED\n");
}

void test_canonicalize_negation() {
    Expr* e = create_neg(create_add(create_var("x"), create_const(2)));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);
    Term* t = c->poly->head;
    ASSERT_COEFF(t, 1, -1); t = t->next;
    ASSERT_COEFF(t, 0, -2); assert(t->next == NULL);
    free_expr(e); free_expr(c);
    printf("  test_canonicalize_negation PASSED\n");
}

void test_canonicalize_power() {
    Expr* e = create_pow(create_add(create_var("x"), create_const(1)), create_const(2));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);
    Term* t = c->poly->head;
    ASSERT_COEFF(t, 2, 1); t = t->next;
    ASSERT_COEFF(t, 1, 2); t = t->next;
    ASSERT_COEFF(t, 0, 1);
    free_expr(e); free_expr(c);
    printf("  test_canonicalize_power PASSED\n");
}

void test_canonicalize_division() {
    Expr* num = create_sub(create_mul(create_var("x"), create_var("x")), create_const(1));
    Expr* den = create_sub(create_var("x"), create_const(1));
    Expr* e = create_div(num, den);
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);
    Term* t = c->poly->head;
    ASSERT_COEFF(t, 1, 1); t = t->next;
    ASSERT_COEFF(t, 0, 1);
    free_expr(e); free_expr(c);
    printf("  test_canonicalize_division PASSED\n");
}

void test_simplify_combines_like_terms() {
    Expr* e = create_add(
        create_mul(create_const(5), create_pow(create_var("x"), create_const(3))),
        create_mul(create_const(3), create_pow(create_var("x"), create_const(3))));
    Expr* s = simplify(e);
    assert(s && s->type == EXPR_POLY);
    ASSERT_COEFF(s->poly->head, 3, 8); assert(s->poly->head->next == NULL);
    free_expr(e); free_expr(s);
    printf("  test_simplify_combines_like_terms PASSED\n");
}

/***********************************************************************
*  factoring tests
*************************************************************************/

void test_factor_linear() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(1); add_term(p, &c, 2); rat_clear(&c);
    c = R(-5); add_term(p, &c, 1); rat_clear(&c);
    c = R(6); add_term(p, &c, 0); rat_clear(&c);
    Factorization* f = factorize(p);
    assert(f && f->count == 2);
    assert(rat_is_one(&f->content));
    Polynomial* product = poly_multiply(f->factors[0], f->factors[1]);
    Term* t = product->head;
    ASSERT_COEFF(t, 2, 1); t = t->next;
    ASSERT_COEFF(t, 1, -5); t = t->next;
    ASSERT_COEFF(t, 0, 6);
    free_polynomial(product); free_factorization(f); free_polynomial(p);
    printf("  test_factor_linear PASSED\n");
}

void test_factor_repeated() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(1); add_term(p, &c, 3); rat_clear(&c);
    c = R(-3); add_term(p, &c, 2); rat_clear(&c);
    c = R(3); add_term(p, &c, 1); rat_clear(&c);
    c = R(-1); add_term(p, &c, 0); rat_clear(&c);
    Factorization* f = factorize(p);
    assert(f && f->count == 1 && f->multiplicities[0] == 3);
    ASSERT_COEFF(f->factors[0]->head, 1, 1);
    ASSERT_COEFF(f->factors[0]->head->next, 0, -1);
    free_factorization(f); free_polynomial(p);
    printf("  test_factor_repeated PASSED\n");
}

void test_factor_with_content() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(6); add_term(p, &c, 2); rat_clear(&c);
    c = R(12); add_term(p, &c, 1); rat_clear(&c);
    c = R(6); add_term(p, &c, 0); rat_clear(&c);
    Factorization* f = factorize(p);
    assert(f);
    Rational six = R(6);
    assert(rat_eq(&f->content, &six));
    rat_clear(&six);
    assert(f->count == 1 && f->multiplicities[0] == 2);
    free_factorization(f); free_polynomial(p);
    printf("  test_factor_with_content PASSED\n");
}

void test_factor_irreducible() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(1); add_term(p, &c, 2); rat_clear(&c);
    c = R(1); add_term(p, &c, 0); rat_clear(&c);
    Factorization* f = factorize(p);
    assert(f && f->count == 1 && f->multiplicities[0] == 1);
    assert(polynomial_degree(f->factors[0]) == 2);
    free_factorization(f); free_polynomial(p);
    printf("  test_factor_irreducible PASSED\n");
}

void test_factor_rational_roots() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(2); add_term(p, &c, 3); rat_clear(&c);
    c = R(-3); add_term(p, &c, 2); rat_clear(&c);
    c = R(-8); add_term(p, &c, 1); rat_clear(&c);
    c = R(-3); add_term(p, &c, 0); rat_clear(&c);
    Factorization* f = factorize(p);
    assert(f);
    int total_degree = 0;
    for (int i = 0; i < f->count; i++)
        total_degree += polynomial_degree(f->factors[i]) * f->multiplicities[i];
    assert(total_degree == 3);
    free_factorization(f); free_polynomial(p);
    printf("  test_factor_rational_roots PASSED\n");
}

void test_factor_to_string() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(1); add_term(p, &c, 2); rat_clear(&c);
    c = R(-5); add_term(p, &c, 1); rat_clear(&c);
    c = R(6); add_term(p, &c, 0); rat_clear(&c);
    Factorization* f = factorize(p);
    char* s = factorization_to_string(f);
    assert(strstr(s, "(") != NULL);
    printf("  factored form: %s\n", s);
    free(s); free_factorization(f); free_polynomial(p);
    printf("  test_factor_to_string PASSED\n");
}

/***********************************************************************
*  solver tests
*************************************************************************/

void test_solve_linear() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(2); add_term(p, &c, 1); rat_clear(&c);
    c = R(4); add_term(p, &c, 0); rat_clear(&c);
    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 1);
    Rational neg2 = R(-2);
    assert(rat_eq(&sol->rational_roots[0], &neg2));
    rat_clear(&neg2);
    free_solutions(sol); free_polynomial(p);
    printf("  test_solve_linear PASSED\n");
}

void test_solve_quadratic_rational() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(1); add_term(p, &c, 2); rat_clear(&c);
    c = R(-5); add_term(p, &c, 1); rat_clear(&c);
    c = R(6); add_term(p, &c, 0); rat_clear(&c);
    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 2);
    Rational two = R(2), three = R(3);
    bool has_2 = false, has_3 = false;
    for (int i = 0; i < sol->num_rational; i++) {
        if (rat_eq(&sol->rational_roots[i], &two)) has_2 = true;
        if (rat_eq(&sol->rational_roots[i], &three)) has_3 = true;
    }
    assert(has_2 && has_3);
    rat_clear(&two); rat_clear(&three);
    free_solutions(sol); free_polynomial(p);
    printf("  test_solve_quadratic_rational PASSED\n");
}

void test_solve_quadratic_irrational() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(1); add_term(p, &c, 2); rat_clear(&c);
    c = R(-2); add_term(p, &c, 0); rat_clear(&c);
    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 0 && sol->num_irrational == 2);
    assert(fabs(fabs(sol->irrational_roots[0]) - sqrt(2.0)) < 1e-10);
    free_solutions(sol); free_polynomial(p);
    printf("  test_solve_quadratic_irrational PASSED\n");
}

void test_solve_cubic() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(1); add_term(p, &c, 3); rat_clear(&c);
    c = R(-6); add_term(p, &c, 2); rat_clear(&c);
    c = R(11); add_term(p, &c, 1); rat_clear(&c);
    c = R(-6); add_term(p, &c, 0); rat_clear(&c);
    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 3);
    char* s = solutions_to_string(sol);
    printf("  solved: %s\n", s);
    free(s); free_solutions(sol); free_polynomial(p);
    printf("  test_solve_cubic PASSED\n");
}

void test_solve_no_real_roots() {
    Polynomial* p = create_polynomial();
    Rational c; c = R(1); add_term(p, &c, 2); rat_clear(&c);
    c = R(1); add_term(p, &c, 0); rat_clear(&c);
    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 0 && sol->num_irrational == 0 && sol->num_irreducible == 1);
    free_solutions(sol); free_polynomial(p);
    printf("  test_solve_no_real_roots PASSED\n");
}

/***********************************************************************
*  calculus tests
*************************************************************************/

void test_differentiate_poly() {
    Expr* e = create_add(
        create_add(create_mul(create_const(3), create_pow(create_var("x"), create_const(2))),
                   create_mul(create_const(4), create_var("x"))),
        create_const(2));
    Expr* d = differentiate(e, "x");
    assert(d && d->type == EXPR_POLY);
    Term* t = d->poly->head;
    ASSERT_COEFF(t, 1, 6); t = t->next;
    ASSERT_COEFF(t, 0, 4); assert(t->next == NULL);
    free_expr(e); free_expr(d);
    printf("  test_differentiate_poly PASSED\n");
}

void test_differentiate_via_parser() {
    Expr* parsed = parse_and_simplify("x^3 + 2x^2 + x");
    Expr* d = differentiate(parsed, "x");
    assert(d && d->type == EXPR_POLY);
    Term* t = d->poly->head;
    ASSERT_COEFF(t, 2, 3); t = t->next;
    ASSERT_COEFF(t, 1, 4); t = t->next;
    ASSERT_COEFF(t, 0, 1); assert(t->next == NULL);
    free_expr(parsed); free_expr(d);
    printf("  test_differentiate_via_parser PASSED\n");
}

void test_integrate_poly() {
    Expr* e = create_add(
        create_mul(create_const(3), create_pow(create_var("x"), create_const(2))),
        create_mul(create_const(2), create_var("x")));
    Expr* i = integrate(e, "x");
    assert(i && i->type == EXPR_POLY);
    Term* t = i->poly->head;
    ASSERT_COEFF(t, 3, 1); t = t->next;
    ASSERT_COEFF(t, 2, 1); assert(t->next == NULL);
    free_expr(e); free_expr(i);
    printf("  test_integrate_poly PASSED\n");
}

/***********************************************************************
*  matrix / linear systems tests
*************************************************************************/

void test_matrix_create_free() {
    Matrix* m = matrix_create(2, 3);
    assert(m->rows == 2 && m->cols == 3);
    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 3; c++)
            assert(rat_is_zero(&MAT_AT(m, r, c)));
    Rational seven = R(7);
    rat_set(&MAT_AT(m, 0, 1), &seven);
    assert(rat_eq(&MAT_AT(m, 0, 1), &seven));
    rat_clear(&seven);
    matrix_free(m);
    printf("  test_matrix_create_free PASSED\n");
}

void test_matrix_copy() {
    Matrix* m = matrix_create(2, 2);
    Rational c;
    c = R(1); rat_set(&MAT_AT(m, 0, 0), &c); rat_clear(&c);
    c = R(2); rat_set(&MAT_AT(m, 0, 1), &c); rat_clear(&c);
    c = R(3); rat_set(&MAT_AT(m, 1, 0), &c); rat_clear(&c);
    c = R(4); rat_set(&MAT_AT(m, 1, 1), &c); rat_clear(&c);
    Matrix* cp = matrix_copy(m);
    c = R(99); rat_set(&MAT_AT(m, 0, 0), &c); rat_clear(&c);
    Rational one = R(1);
    assert(rat_eq(&MAT_AT(cp, 0, 0), &one));
    rat_clear(&one);
    matrix_free(m); matrix_free(cp);
    printf("  test_matrix_copy PASSED\n");
}

void test_vector_create_free() {
    Vector* v = vector_create(3);
    Rational c = R(5); rat_set(&v->data[0], &c); rat_clear(&c);
    c = rat_from_si(1, 2); rat_set(&v->data[2], &c);
    assert(rat_eq(&v->data[2], &c));
    rat_clear(&c);
    vector_free(v);
    printf("  test_vector_create_free PASSED\n");
}

void test_matrix_to_string() {
    Matrix* m = matrix_create(2, 2);
    Rational c;
    c = R(1); rat_set(&MAT_AT(m, 0, 0), &c); rat_clear(&c);
    c = R(4); rat_set(&MAT_AT(m, 1, 1), &c); rat_clear(&c);
    char* s = matrix_to_string(m);
    assert(strstr(s, "1") && strstr(s, "4"));
    free(s); matrix_free(m);
    printf("  test_matrix_to_string PASSED\n");
}

void test_rref_identity() {
    Matrix* m = matrix_create(3, 3);
    for (int i = 0; i < 3; i++) { Rational one = R(1); rat_set(&MAT_AT(m, i, i), &one); rat_clear(&one); }
    Matrix* rref = matrix_rref(m);
    Rational one = rat_one_val(), zero = rat_zero_val();
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            assert(rat_eq(&MAT_AT(rref, r, c), (r == c) ? &one : &zero));
    rat_clear(&one); rat_clear(&zero);
    matrix_free(m); matrix_free(rref);
    printf("  test_rref_identity PASSED\n");
}

void test_rref_2x3() {
    Matrix* m = matrix_create(2, 3);
    Rational vals[] = {R(1),R(2),R(3),R(4),R(5),R(6)};
    for (int i = 0; i < 6; i++) { rat_set(&m->data[i], &vals[i]); rat_clear(&vals[i]); }
    Matrix* rref = matrix_rref(m);
    Rational expected[] = {R(1),R(0),R(-1),R(0),R(1),R(2)};
    for (int i = 0; i < 6; i++) { assert(rat_eq(&rref->data[i], &expected[i])); rat_clear(&expected[i]); }
    matrix_free(m); matrix_free(rref);
    printf("  test_rref_2x3 PASSED\n");
}

void test_rref_needs_swap() {
    Matrix* m = matrix_create(2, 2);
    Rational c;
    c = R(0); rat_set(&MAT_AT(m, 0, 0), &c); rat_clear(&c);
    c = R(1); rat_set(&MAT_AT(m, 0, 1), &c); rat_clear(&c);
    c = R(1); rat_set(&MAT_AT(m, 1, 0), &c); rat_clear(&c);
    c = R(0); rat_set(&MAT_AT(m, 1, 1), &c); rat_clear(&c);
    Matrix* rref = matrix_rref(m);
    Rational one = R(1);
    assert(rat_eq(&MAT_AT(rref, 0, 0), &one));
    assert(rat_eq(&MAT_AT(rref, 1, 1), &one));
    rat_clear(&one);
    matrix_free(m); matrix_free(rref);
    printf("  test_rref_needs_swap PASSED\n");
}

void test_rank() {
    Matrix* m = matrix_create(3, 3);
    Rational vals[] = {R(1),R(2),R(3),R(4),R(5),R(6),R(7),R(8),R(9)};
    for (int i = 0; i < 9; i++) { rat_set(&m->data[i], &vals[i]); rat_clear(&vals[i]); }
    assert(matrix_rank(m) == 2);
    matrix_free(m);
    printf("  test_rank PASSED\n");
}

void test_det_2x2() {
    Matrix* m = matrix_create(2, 2);
    Rational vals[] = {R(1),R(2),R(3),R(4)};
    for (int i = 0; i < 4; i++) { rat_set(&m->data[i], &vals[i]); rat_clear(&vals[i]); }
    Rational d; rat_init(&d);
    matrix_det(&d, m);
    Rational neg2 = R(-2);
    assert(rat_eq(&d, &neg2));
    rat_clear(&neg2); rat_clear(&d);
    matrix_free(m);
    printf("  test_det_2x2 PASSED\n");
}

void test_det_3x3() {
    Matrix* m = matrix_create(3, 3);
    Rational vals[] = {R(2),R(1),R(1),R(1),R(3),R(2),R(1),R(0),R(0)};
    for (int i = 0; i < 9; i++) { rat_set(&m->data[i], &vals[i]); rat_clear(&vals[i]); }
    Rational d; rat_init(&d);
    matrix_det(&d, m);
    Rational neg1 = R(-1);
    assert(rat_eq(&d, &neg1));
    rat_clear(&neg1); rat_clear(&d);
    matrix_free(m);
    printf("  test_det_3x3 PASSED\n");
}

void test_det_singular() {
    Matrix* m = matrix_create(2, 2);
    Rational vals[] = {R(1),R(2),R(2),R(4)};
    for (int i = 0; i < 4; i++) { rat_set(&m->data[i], &vals[i]); rat_clear(&vals[i]); }
    Rational d; rat_init(&d);
    matrix_det(&d, m);
    assert(rat_is_zero(&d));
    rat_clear(&d); matrix_free(m);
    printf("  test_det_singular PASSED\n");
}

void test_det_with_swap() {
    Matrix* m = matrix_create(2, 2);
    Rational vals[] = {R(0),R(1),R(1),R(0)};
    for (int i = 0; i < 4; i++) { rat_set(&m->data[i], &vals[i]); rat_clear(&vals[i]); }
    Rational d; rat_init(&d);
    matrix_det(&d, m);
    Rational neg1 = R(-1);
    assert(rat_eq(&d, &neg1));
    rat_clear(&neg1); rat_clear(&d); matrix_free(m);
    printf("  test_det_with_swap PASSED\n");
}

void test_solve_2x2_unique() {
    Matrix* A = matrix_create(2, 2);
    Rational c;
    c = R(1); rat_set(&MAT_AT(A, 0, 0), &c); rat_clear(&c);
    c = R(2); rat_set(&MAT_AT(A, 0, 1), &c); rat_clear(&c);
    c = R(3); rat_set(&MAT_AT(A, 1, 0), &c); rat_clear(&c);
    c = R(4); rat_set(&MAT_AT(A, 1, 1), &c); rat_clear(&c);
    Vector* b = vector_create(2);
    c = R(5); rat_set(&b->data[0], &c); rat_clear(&c);
    c = R(11); rat_set(&b->data[1], &c); rat_clear(&c);
    LinSysResult* res = linsys_solve(A, b);
    assert(res->status == LINSYS_UNIQUE);
    Rational e1 = R(1), e2 = R(2);
    assert(rat_eq(&res->solution->data[0], &e1));
    assert(rat_eq(&res->solution->data[1], &e2));
    rat_clear(&e1); rat_clear(&e2);
    linsys_result_free(res); matrix_free(A); vector_free(b);
    printf("  test_solve_2x2_unique PASSED\n");
}

void test_solve_3x3_unique() {
    Matrix* A = matrix_create(3, 3);
    Rational c;
    c = R(1); rat_set(&MAT_AT(A, 0, 0), &c); rat_clear(&c);
    c = R(1); rat_set(&MAT_AT(A, 0, 1), &c); rat_clear(&c);
    c = R(1); rat_set(&MAT_AT(A, 0, 2), &c); rat_clear(&c);
    c = R(0); rat_set(&MAT_AT(A, 1, 0), &c); rat_clear(&c);
    c = R(2); rat_set(&MAT_AT(A, 1, 1), &c); rat_clear(&c);
    c = R(5); rat_set(&MAT_AT(A, 1, 2), &c); rat_clear(&c);
    c = R(2); rat_set(&MAT_AT(A, 2, 0), &c); rat_clear(&c);
    c = R(5); rat_set(&MAT_AT(A, 2, 1), &c); rat_clear(&c);
    c = R(-1); rat_set(&MAT_AT(A, 2, 2), &c); rat_clear(&c);
    Vector* b = vector_create(3);
    c = R(6); rat_set(&b->data[0], &c); rat_clear(&c);
    c = R(-4); rat_set(&b->data[1], &c); rat_clear(&c);
    c = R(27); rat_set(&b->data[2], &c); rat_clear(&c);
    LinSysResult* res = linsys_solve(A, b);
    assert(res->status == LINSYS_UNIQUE);
    Rational e1 = R(5), e2 = R(3), e3 = R(-2);
    assert(rat_eq(&res->solution->data[0], &e1));
    assert(rat_eq(&res->solution->data[1], &e2));
    assert(rat_eq(&res->solution->data[2], &e3));
    rat_clear(&e1); rat_clear(&e2); rat_clear(&e3);
    linsys_result_free(res); matrix_free(A); vector_free(b);
    printf("  test_solve_3x3_unique PASSED\n");
}

void test_solve_inconsistent() {
    Matrix* A = matrix_create(2, 2);
    Rational c;
    c = R(1); rat_set(&MAT_AT(A, 0, 0), &c); rat_clear(&c);
    c = R(1); rat_set(&MAT_AT(A, 0, 1), &c); rat_clear(&c);
    c = R(1); rat_set(&MAT_AT(A, 1, 0), &c); rat_clear(&c);
    c = R(1); rat_set(&MAT_AT(A, 1, 1), &c); rat_clear(&c);
    Vector* b = vector_create(2);
    c = R(1); rat_set(&b->data[0], &c); rat_clear(&c);
    c = R(2); rat_set(&b->data[1], &c); rat_clear(&c);
    LinSysResult* res = linsys_solve(A, b);
    assert(res->status == LINSYS_INCONSISTENT);
    linsys_result_free(res); matrix_free(A); vector_free(b);
    printf("  test_solve_inconsistent PASSED\n");
}

void test_solve_infinite() {
    Matrix* A = matrix_create(1, 3);
    Rational c;
    c = R(1); rat_set(&MAT_AT(A, 0, 0), &c); rat_clear(&c);
    c = R(2); rat_set(&MAT_AT(A, 0, 1), &c); rat_clear(&c);
    c = R(3); rat_set(&MAT_AT(A, 0, 2), &c); rat_clear(&c);
    Vector* b = vector_create(1);
    c = R(1); rat_set(&b->data[0], &c); rat_clear(&c);
    LinSysResult* res = linsys_solve(A, b);
    assert(res->status == LINSYS_INFINITE && res->rank == 1);
    linsys_result_free(res); matrix_free(A); vector_free(b);
    printf("  test_solve_infinite PASSED\n");
}

void test_solve_hilbert_3x3() {
    Matrix* H = matrix_create(3, 3);
    Rational c;
    c = rat_from_si(1,1); rat_set(&MAT_AT(H,0,0), &c); rat_clear(&c);
    c = rat_from_si(1,2); rat_set(&MAT_AT(H,0,1), &c); rat_clear(&c);
    c = rat_from_si(1,3); rat_set(&MAT_AT(H,0,2), &c); rat_clear(&c);
    c = rat_from_si(1,2); rat_set(&MAT_AT(H,1,0), &c); rat_clear(&c);
    c = rat_from_si(1,3); rat_set(&MAT_AT(H,1,1), &c); rat_clear(&c);
    c = rat_from_si(1,4); rat_set(&MAT_AT(H,1,2), &c); rat_clear(&c);
    c = rat_from_si(1,3); rat_set(&MAT_AT(H,2,0), &c); rat_clear(&c);
    c = rat_from_si(1,4); rat_set(&MAT_AT(H,2,1), &c); rat_clear(&c);
    c = rat_from_si(1,5); rat_set(&MAT_AT(H,2,2), &c); rat_clear(&c);
    Vector* b = vector_create(3);
    Rational one = rat_one_val();
    rat_set(&b->data[0], &one); rat_set(&b->data[1], &one); rat_set(&b->data[2], &one);
    rat_clear(&one);
    LinSysResult* res = linsys_solve(H, b);
    assert(res->status == LINSYS_UNIQUE);
    Rational e1 = R(3), e2 = R(-24), e3 = R(30);
    assert(rat_eq(&res->solution->data[0], &e1));
    assert(rat_eq(&res->solution->data[1], &e2));
    assert(rat_eq(&res->solution->data[2], &e3));
    rat_clear(&e1); rat_clear(&e2); rat_clear(&e3);
    char* s = linsys_result_to_string(res);
    printf("  hilbert solved: %s\n", s);
    free(s);
    linsys_result_free(res); matrix_free(H); vector_free(b);
    printf("  test_solve_hilbert_3x3 PASSED\n");
}

void test_inverse_2x2() {
    Matrix* A = matrix_create(2, 2);
    Rational vals[] = {R(1),R(2),R(3),R(4)};
    for (int i = 0; i < 4; i++) { rat_set(&A->data[i], &vals[i]); rat_clear(&vals[i]); }
    Matrix* inv = matrix_inverse(A);
    assert(inv != NULL);
    Rational expected[] = {R(-2),R(1),rat_from_si(3,2),rat_from_si(-1,2)};
    for (int i = 0; i < 4; i++) { assert(rat_eq(&inv->data[i], &expected[i])); rat_clear(&expected[i]); }
    matrix_free(A); matrix_free(inv);
    printf("  test_inverse_2x2 PASSED\n");
}

void test_inverse_singular() {
    Matrix* m = matrix_create(2, 2);
    Rational vals[] = {R(1),R(2),R(2),R(4)};
    for (int i = 0; i < 4; i++) { rat_set(&m->data[i], &vals[i]); rat_clear(&vals[i]); }
    assert(matrix_inverse(m) == NULL);
    matrix_free(m);
    printf("  test_inverse_singular PASSED\n");
}

void test_matrix_parse() {
    Matrix* m = matrix_parse("[[1, 2], [3, 4]]");
    assert(m && m->rows == 2 && m->cols == 2);
    Rational one = R(1), four = R(4);
    assert(rat_eq(&MAT_AT(m, 0, 0), &one));
    assert(rat_eq(&MAT_AT(m, 1, 1), &four));
    rat_clear(&one); rat_clear(&four);
    matrix_free(m);

    Matrix* m2 = matrix_parse("[[1/2, -3/4], [0, 5]]");
    assert(m2 && m2->rows == 2 && m2->cols == 2);
    Rational half = rat_from_si(1, 2);
    assert(rat_eq(&MAT_AT(m2, 0, 0), &half));
    rat_clear(&half);
    matrix_free(m2);

    printf("  test_matrix_parse PASSED\n");
}

void test_vector_parse() {
    Vector* v = vector_parse("[1, -2, 3/5]");
    assert(v && v->n == 3);
    Rational e1 = R(1), e2 = R(-2), e3 = rat_from_si(3, 5);
    assert(rat_eq(&v->data[0], &e1));
    assert(rat_eq(&v->data[1], &e2));
    assert(rat_eq(&v->data[2], &e3));
    rat_clear(&e1); rat_clear(&e2); rat_clear(&e3);
    vector_free(v);
    assert(vector_parse("not a vector") == NULL);
    printf("  test_vector_parse PASSED\n");
}

/***********************************************************************
*  partial fraction tests
*************************************************************************/

void test_pf_simple() {
    Polynomial* num = create_polynomial();
    Rational c; c = R(2); add_term(num, &c, 1); rat_clear(&c);
    c = R(3); add_term(num, &c, 0); rat_clear(&c);
    Polynomial* den = create_polynomial();
    c = R(1); add_term(den, &c, 2); rat_clear(&c);
    c = R(-1); add_term(den, &c, 0); rat_clear(&c);
    PartialFractionDecomp* pf = partial_fractions(num, den);
    assert(pf && pf->num_terms == 2);
    char* s = pf_to_string(pf);
    printf("  (2x+3)/(x^2-1) = %s\n", s);
    free(s);
    // verify coefficients
    Rational half52 = rat_from_si(5, 2), neg_half = rat_from_si(-1, 2);
    for (int i = 0; i < pf->num_terms; i++) {
        Rational root_neg;
        rat_init(&root_neg);
        test_get_coeff(&root_neg, pf->terms[i].factor, 0);
        Rational neg1 = R(-1);
        if (rat_eq(&root_neg, &neg1))
            assert(rat_eq(&pf->terms[i].coeffs[0], &half52));
        else
            assert(rat_eq(&pf->terms[i].coeffs[0], &neg_half));
        rat_clear(&neg1); rat_clear(&root_neg);
    }
    rat_clear(&half52); rat_clear(&neg_half);
    free_pf_decomp(pf); free_polynomial(num); free_polynomial(den);
    printf("  test_pf_simple PASSED\n");
}

void test_pf_repeated_root() {
    Polynomial* num = create_polynomial();
    Rational c; c = R(3); add_term(num, &c, 1); rat_clear(&c);
    c = R(5); add_term(num, &c, 0); rat_clear(&c);
    Polynomial* den = create_polynomial();
    c = R(1); add_term(den, &c, 2); rat_clear(&c);
    c = R(-2); add_term(den, &c, 1); rat_clear(&c);
    c = R(1); add_term(den, &c, 0); rat_clear(&c);
    PartialFractionDecomp* pf = partial_fractions(num, den);
    assert(pf && pf->num_terms == 1 && pf->terms[0].power == 2);
    Rational e3 = R(3), e8 = R(8);
    assert(rat_eq(&pf->terms[0].coeffs[0], &e3));
    assert(rat_eq(&pf->terms[0].coeffs[1], &e8));
    rat_clear(&e3); rat_clear(&e8);
    char* s = pf_to_string(pf);
    printf("  (3x+5)/(x-1)^2 = %s\n", s);
    free(s);
    free_pf_decomp(pf); free_polynomial(num); free_polynomial(den);
    printf("  test_pf_repeated_root PASSED\n");
}

void test_pf_improper() {
    Polynomial* num = create_polynomial();
    Rational c; c = R(1); add_term(num, &c, 2); rat_clear(&c);
    c = R(1); add_term(num, &c, 0); rat_clear(&c);
    Polynomial* den = create_polynomial();
    c = R(1); add_term(den, &c, 2); rat_clear(&c);
    c = R(-1); add_term(den, &c, 0); rat_clear(&c);
    PartialFractionDecomp* pf = partial_fractions(num, den);
    assert(pf && pf->poly_part && !poly_is_zero(pf->poly_part));
    ASSERT_COEFF(pf->poly_part->head, 0, 1);
    assert(pf->num_terms == 2);
    char* s = pf_to_string(pf);
    printf("  (x^2+1)/(x^2-1) = %s\n", s);
    free(s);
    free_pf_decomp(pf); free_polynomial(num); free_polynomial(den);
    printf("  test_pf_improper PASSED\n");
}

/***********************************************************************
*  rational function integration tests
*************************************************************************/

void test_integrate_rational_simple() {
    Polynomial* num = create_polynomial();
    Rational c; c = R(2); add_term(num, &c, 1); rat_clear(&c);
    c = R(3); add_term(num, &c, 0); rat_clear(&c);
    Polynomial* den = create_polynomial();
    c = R(1); add_term(den, &c, 2); rat_clear(&c);
    c = R(-1); add_term(den, &c, 0); rat_clear(&c);
    Expr* ratfn = create_rational_fn(num, den);
    Expr* result = integrate(ratfn, "x");
    assert(result != NULL);
    char* s = expr_to_string(result);
    printf("  int(2x+3)/(x²-1)dx = %s\n", s);
    assert(strstr(s, "ln") != NULL);
    free(s); free_expr(ratfn); free_expr(result);
    printf("  test_integrate_rational_simple PASSED\n");
}

void test_integrate_rational_repeated() {
    Polynomial* num = create_polynomial();
    Rational c; c = R(1); add_term(num, &c, 0); rat_clear(&c);
    Polynomial* den = create_polynomial();
    c = R(1); add_term(den, &c, 2); rat_clear(&c);
    c = R(-2); add_term(den, &c, 1); rat_clear(&c);
    c = R(1); add_term(den, &c, 0); rat_clear(&c);
    Expr* ratfn = create_rational_fn(num, den);
    Expr* result = integrate(ratfn, "x");
    assert(result != NULL);
    char* s = expr_to_string(result);
    printf("  int 1/(x-1)²dx = %s\n", s);
    assert(strstr(s, "ln") == NULL);
    free(s); free_expr(ratfn); free_expr(result);
    printf("  test_integrate_rational_repeated PASSED\n");
}

void test_integrate_rational_improper() {
    Polynomial* num = create_polynomial();
    Rational c; c = R(1); add_term(num, &c, 2); rat_clear(&c);
    c = R(1); add_term(num, &c, 0); rat_clear(&c);
    Polynomial* den = create_polynomial();
    c = R(1); add_term(den, &c, 2); rat_clear(&c);
    c = R(-1); add_term(den, &c, 0); rat_clear(&c);
    Expr* ratfn = create_rational_fn(num, den);
    Expr* result = integrate(ratfn, "x");
    assert(result != NULL);
    char* s = expr_to_string(result);
    printf("  int(x²+1)/(x²-1)dx = %s\n", s);
    assert(strstr(s, "ln") != NULL);
    free(s); free_expr(ratfn); free_expr(result);
    printf("  test_integrate_rational_improper PASSED\n");
}

/***********************************************************************
*  complex / stress tests (GMP-ported)
*************************************************************************/

void test_parse_expand_product() {
    Expr* e = parse_and_simplify("(x + 1) * (x + 2) * (x + 3)");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    Rational e1 = R(1); assert(t && t->exponent == 3 && rat_eq(&t->coeff, &e1)); rat_clear(&e1);
    t = t->next;
    Rational e6 = R(6); assert(t && t->exponent == 2 && rat_eq(&t->coeff, &e6)); rat_clear(&e6);
    t = t->next;
    Rational e11 = R(11); assert(t && t->exponent == 1 && rat_eq(&t->coeff, &e11)); rat_clear(&e11);
    t = t->next;
    e6 = R(6); assert(t && t->exponent == 0 && rat_eq(&t->coeff, &e6)); rat_clear(&e6);
    assert(t->next == NULL);
    free_expr(e);
    printf("  test_parse_expand_product PASSED\n");
}

void test_expand_then_factor_roundtrip() {
    // (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
    Polynomial* p = create_polynomial();
    Rational c;
    c = R(1); add_term(p, &c, 4); rat_clear(&c);
    c = R(-10); add_term(p, &c, 3); rat_clear(&c);
    c = R(35); add_term(p, &c, 2); rat_clear(&c);
    c = R(-50); add_term(p, &c, 1); rat_clear(&c);
    c = R(24); add_term(p, &c, 0); rat_clear(&c);
    Factorization* f = factorize(p);
    assert(f && f->count == 4);
    for (int i = 0; i < 4; i++) {
        assert(f->multiplicities[i] == 1);
        assert(polynomial_degree(f->factors[i]) == 1);
    }
    free_factorization(f);
    free_polynomial(p);
    printf("  test_expand_then_factor_roundtrip PASSED\n");
}

void test_differentiate_integrate_inverse() {
    Polynomial* fp = create_polynomial();
    Rational c;
    c = R(4); add_term(fp, &c, 3); rat_clear(&c);
    c = R(-6); add_term(fp, &c, 2); rat_clear(&c);
    c = R(2); add_term(fp, &c, 1); rat_clear(&c);
    c = R(-7); add_term(fp, &c, 0); rat_clear(&c);
    Polynomial* integral = poly_integral(fp);
    Polynomial* deriv = poly_derivative(integral);
    Term* tf = fp->head; Term* td = deriv->head;
    while (tf && td) {
        assert(tf->exponent == td->exponent);
        assert(rat_eq(&tf->coeff, &td->coeff));
        tf = tf->next; td = td->next;
    }
    assert(tf == NULL && td == NULL);
    free_polynomial(fp); free_polynomial(integral); free_polynomial(deriv);
    printf("  test_differentiate_integrate_inverse PASSED\n");
}

void test_high_degree_polynomial_multiply() {
    // (x^5 + 1)(x^5 - 1) = x^10 - 1
    Polynomial* a = create_polynomial();
    Polynomial* b = create_polynomial();
    Rational one = R(1), neg1 = R(-1);
    add_term(a, &one, 5); add_term(a, &one, 0);
    add_term(b, &one, 5); add_term(b, &neg1, 0);
    Polynomial* product = poly_multiply(a, b);
    Term* t = product->head;
    assert(t && t->exponent == 10 && rat_eq(&t->coeff, &one));
    t = t->next;
    assert(t && t->exponent == 0 && rat_eq(&t->coeff, &neg1));
    assert(t->next == NULL);
    rat_clear(&one); rat_clear(&neg1);
    free_polynomial(a); free_polynomial(b); free_polynomial(product);
    printf("  test_high_degree_polynomial_multiply PASSED\n");
}

void test_solve_quartic() {
    // x^4 - 5x^2 + 4 = (x-1)(x+1)(x-2)(x+2)
    Polynomial* p = create_polynomial();
    Rational c;
    c = R(1); add_term(p, &c, 4); rat_clear(&c);
    c = R(-5); add_term(p, &c, 2); rat_clear(&c);
    c = R(4); add_term(p, &c, 0); rat_clear(&c);
    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 4);
    free_solutions(sol); free_polynomial(p);
    printf("  test_solve_quartic PASSED\n");
}

void test_solve_degree5_partial() {
    // x^5 - x = x(x-1)(x+1)(x^2+1)
    Polynomial* p = create_polynomial();
    Rational c;
    c = R(1); add_term(p, &c, 5); rat_clear(&c);
    c = R(-1); add_term(p, &c, 1); rat_clear(&c);
    Solutions* sol = solve_polynomial(p);
    assert(sol->num_rational == 3);
    assert(sol->num_irreducible == 1);
    free_solutions(sol); free_polynomial(p);
    printf("  test_solve_degree5_partial PASSED\n");
}

void test_parse_complex_expression() {
    Expr* e = parse_and_simplify("(2x^2 + 3x - 5) * (x - 1) - (x^3 + 2x^2 - 8x + 5)");
    assert(e && e->type == EXPR_POLY);
    Term* t = e->poly->head;
    Rational e1 = R(1), en1 = R(-1);
    assert(t && t->exponent == 3 && rat_eq(&t->coeff, &e1));
    t = t->next;
    assert(t && t->exponent == 2 && rat_eq(&t->coeff, &en1));
    assert(t->next == NULL);
    rat_clear(&e1); rat_clear(&en1);
    free_expr(e);
    printf("  test_parse_complex_expression PASSED\n");
}

void test_rational_arithmetic_stress() {
    // ((1/7 + 1/11) * 77) should be exactly 18
    Rational a = rat_from_si(1, 7);
    Rational b = rat_from_si(1, 11);
    Rational sum; rat_init(&sum);
    rat_add(&sum, &a, &b);
    Rational seventy7 = R(77);
    Rational result; rat_init(&result);
    rat_mul(&result, &sum, &seventy7);
    Rational e18 = R(18);
    assert(rat_eq(&result, &e18));
    rat_clear(&a); rat_clear(&b); rat_clear(&sum);
    rat_clear(&seventy7); rat_clear(&result); rat_clear(&e18);

    // (1/3)^10 * 3^10 = 1
    Rational third = rat_from_si(1, 3);
    Rational t10; rat_init(&t10);
    rat_pow_int(&t10, &third, 10);
    Rational three = R(3);
    Rational three10; rat_init(&three10);
    rat_pow_int(&three10, &three, 10);
    Rational product; rat_init(&product);
    rat_mul(&product, &t10, &three10);
    Rational one = rat_one_val();
    assert(rat_eq(&product, &one));
    rat_clear(&third); rat_clear(&t10); rat_clear(&three);
    rat_clear(&three10); rat_clear(&product); rat_clear(&one);

    printf("  test_rational_arithmetic_stress PASSED\n");
}

void test_power_of_binomial() {
    // (x + 2)^5 = x^5 + 10x^4 + 40x^3 + 80x^2 + 80x + 32
    Expr* e = create_pow(
        create_add(create_var("x"), create_const(2)),
        create_const(5));
    Expr* c = canonicalize(e);
    assert(c && c->type == EXPR_POLY);
    Term* t = c->poly->head;
    Rational expected;
    expected = R(1); assert(t && t->exponent == 5 && rat_eq(&t->coeff, &expected)); rat_clear(&expected);
    t = t->next;
    expected = R(10); assert(t && t->exponent == 4 && rat_eq(&t->coeff, &expected)); rat_clear(&expected);
    t = t->next;
    expected = R(40); assert(t && t->exponent == 3 && rat_eq(&t->coeff, &expected)); rat_clear(&expected);
    t = t->next;
    expected = R(80); assert(t && t->exponent == 2 && rat_eq(&t->coeff, &expected)); rat_clear(&expected);
    t = t->next;
    expected = R(80); assert(t && t->exponent == 1 && rat_eq(&t->coeff, &expected)); rat_clear(&expected);
    t = t->next;
    expected = R(32); assert(t && t->exponent == 0 && rat_eq(&t->coeff, &expected)); rat_clear(&expected);
    assert(t->next == NULL);
    free_expr(e); free_expr(c);
    printf("  test_power_of_binomial PASSED\n");
}

void test_gmp_large_coefficients() {
    // (1000000007)^3 exceeds 2^63
    Rational big = rat_from_int(1000000007);
    Rational big_cubed; rat_init(&big_cubed);
    rat_pow_int(&big_cubed, &big, 3);
    assert(!rat_is_zero(&big_cubed));
    char* s = rat_to_string(&big_cubed);
    assert(strlen(s) > 18);
    printf("  big^3 = %s\n", s);
    free(s);
    rat_clear(&big); rat_clear(&big_cubed);
    printf("  test_gmp_large_coefficients PASSED\n");
}

/***********************************************************************
*  expression formatter tests
*************************************************************************/

static void assert_expr_str(Expr* e, const char* expected) {
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: expected \"%s\", got \"%s\"\n", expected, s);
        assert(0);
    }
    free(s);
    free_expr(e);
}

void test_fmt_add_negative_rhs() {
    // 3 + (-5) → "3 - 5"
    assert_expr_str(
        create_add(create_const(3), create_const(-5)),
        "3 - 5");
    printf("  test_fmt_add_negative_rhs PASSED\n");
}

void test_fmt_mul_coeff_one() {
    // 1 * ln|x| → "ln|x|"
    assert_expr_str(
        create_mul(create_const(1), create_ln(create_var("x"))),
        "ln|x|");
    printf("  test_fmt_mul_coeff_one PASSED\n");
}

void test_fmt_mul_coeff_neg_one() {
    // -1 * ln|x| → "-ln|x|"
    assert_expr_str(
        create_mul(create_const(-1), create_ln(create_var("x"))),
        "-ln|x|");
    printf("  test_fmt_mul_coeff_neg_one PASSED\n");
}

void test_fmt_neg_coeff_in_sum() {
    // a + (-1/2 * b) → "a - 1/2 * b"
    Rational half = rat_from_si(-1, 2);
    assert_expr_str(
        create_add(create_var("a"), create_mul(create_const_rat(&half), create_var("b"))),
        "a - 1/2 * b");
    rat_clear(&half);
    printf("  test_fmt_neg_coeff_in_sum PASSED\n");
}

void test_fmt_pow_neg_one() {
    // x^(-1) → "1/(x)"
    assert_expr_str(
        create_pow(create_var("x"), create_const(-1)),
        "1/(x)");
    printf("  test_fmt_pow_neg_one PASSED\n");
}

void test_fmt_pow_neg_two() {
    // x^(-2) → "1/(x)^2"
    assert_expr_str(
        create_pow(create_var("x"), create_const(-2)),
        "1/(x)^2");
    printf("  test_fmt_pow_neg_two PASSED\n");
}

void test_fmt_div_parens() {
    // (a + b) / (c * d) → "(a + b) / (c * d)"
    assert_expr_str(
        create_div(
            create_add(create_var("a"), create_var("b")),
            create_mul(create_var("c"), create_var("d"))),
        "(a + b) / (c * d)");
    printf("  test_fmt_div_parens PASSED\n");
}

void test_fmt_pow_poly_base() {
    // (x + 1)^2 — polynomial base gets parens
    Polynomial* p = create_polynomial();
    Rational one = rat_one_val();
    add_term(p, &one, 1);
    add_term(p, &one, 0);
    rat_clear(&one);
    assert_expr_str(
        create_pow(create_poly_expr(p), create_const(2)),
        "(x + 1)^2");
    printf("  test_fmt_pow_poly_base PASSED\n");
}

void test_fmt_chained_add_negatives() {
    // a + (-b) + (-c) → "a - b - c"
    assert_expr_str(
        create_add(
            create_add(create_var("a"), create_neg(create_var("b"))),
            create_neg(create_var("c"))),
        "a - b - c");
    printf("  test_fmt_chained_add_negatives PASSED\n");
}

void test_fmt_integration_output() {
    // build the tree that integrate produces for (2x+3)/(x^2-1):
    // 5/2 * ln|x - 1| + (-1/2) * ln|x + 1|
    // should render as: "5/2 * ln|x - 1| - 1/2 * ln|x + 1|"
    Polynomial* p1 = create_polynomial();
    Rational one = rat_one_val(), neg1r = rat_from_int(-1);
    add_term(p1, &one, 1);
    add_term(p1, &neg1r, 0);
    Polynomial* p2 = create_polynomial();
    Rational one2 = rat_one_val();
    add_term(p2, &one2, 1);
    add_term(p2, &one2, 0);

    Rational coeff1 = rat_from_si(5, 2);
    Rational coeff2 = rat_from_si(-1, 2);
    Expr* term1 = create_mul(create_const_rat(&coeff1), create_ln(create_poly_expr(p1)));
    Expr* term2 = create_mul(create_const_rat(&coeff2), create_ln(create_poly_expr(p2)));
    assert_expr_str(
        create_add(term1, term2),
        "5/2 * ln|x - 1| - 1/2 * ln|x + 1|");
    rat_clear(&one); rat_clear(&neg1r); rat_clear(&one2);
    rat_clear(&coeff1); rat_clear(&coeff2);
    printf("  test_fmt_integration_output PASSED\n");
}

/***********************************************************************
*  GMP precision tests — cases that fail with long long or doubles
*************************************************************************/

void test_gmp_factorial_fraction() {
    // 1/20! is an extremely small fraction; verify it survives a roundtrip
    // 20! = 2432902008176640000 (fits in long long but barely)
    // 25! = 15511210043330985984000000 (does NOT fit in long long)
    // compute 1/25! via repeated division
    Rational r = rat_one_val();
    for (int i = 1; i <= 25; i++) {
        Rational divisor = rat_from_int(i);
        Rational tmp; rat_init(&tmp);
        rat_div(&tmp, &r, &divisor);
        rat_set(&r, &tmp);
        rat_clear(&tmp); rat_clear(&divisor);
    }
    // now multiply by 25! to get back to 1
    for (int i = 1; i <= 25; i++) {
        Rational factor = rat_from_int(i);
        Rational tmp; rat_init(&tmp);
        rat_mul(&tmp, &r, &factor);
        rat_set(&r, &tmp);
        rat_clear(&tmp); rat_clear(&factor);
    }
    Rational one = rat_one_val();
    assert(rat_eq(&r, &one));
    rat_clear(&r); rat_clear(&one);
    printf("  test_gmp_factorial_fraction PASSED\n");
}

void test_gmp_hilbert_5x5() {
    // 5x5 Hilbert matrix has condition number ~943656
    // floating-point solvers produce garbage; exact arithmetic is exact
    // H[i][j] = 1/(i+j+1), b = [1,1,1,1,1]
    // exact solution: [5, -120, 630, -1120, 630]
    Matrix* H = matrix_create(5, 5);
    for (int r = 0; r < 5; r++)
        for (int c = 0; c < 5; c++) {
            Rational v = rat_from_si(1, r + c + 1);
            rat_set(&MAT_AT(H, r, c), &v);
            rat_clear(&v);
        }
    Vector* b = vector_create(5);
    for (int i = 0; i < 5; i++) {
        Rational one = rat_one_val();
        rat_set(&b->data[i], &one);
        rat_clear(&one);
    }

    LinSysResult* res = linsys_solve(H, b);
    assert(res->status == LINSYS_UNIQUE);

    long expected[] = {5, -120, 630, -1120, 630};
    for (int i = 0; i < 5; i++) {
        Rational e = rat_from_int(expected[i]);
        assert(rat_eq(&res->solution->data[i], &e));
        rat_clear(&e);
    }

    char* s = linsys_result_to_string(res);
    printf("  hilbert 5x5: %s\n", s);
    free(s);

    linsys_result_free(res); matrix_free(H); vector_free(b);
    printf("  test_gmp_hilbert_5x5 PASSED\n");
}

void test_gmp_poly_gcd_large() {
    // GCD of two degree-10 polynomials with large coefficients
    // p = (x+1000000007)^5, q = (x+1000000007)^3 * (x+13)^2
    // gcd should be (x+1000000007)^3
    Polynomial* factor1 = create_polynomial();
    Rational one = rat_one_val();
    Rational big = rat_from_int(1000000007);
    add_term(factor1, &one, 1);
    add_term(factor1, &big, 0);

    Polynomial* factor2 = create_polynomial();
    Rational thirteen = rat_from_int(13);
    add_term(factor2, &one, 1);
    add_term(factor2, &thirteen, 0);

    // build (x+1000000007)^5
    Polynomial* p = poly_copy(factor1);
    for (int i = 0; i < 4; i++) {
        Polynomial* tmp = poly_multiply(p, factor1);
        free_polynomial(p); p = tmp;
    }

    // build (x+1000000007)^3 * (x+13)^2
    Polynomial* q = poly_copy(factor1);
    for (int i = 0; i < 2; i++) {
        Polynomial* tmp = poly_multiply(q, factor1);
        free_polynomial(q); q = tmp;
    }
    for (int i = 0; i < 2; i++) {
        Polynomial* tmp = poly_multiply(q, factor2);
        free_polynomial(q); q = tmp;
    }

    Polynomial* g = poly_gcd(p, q);
    // gcd should be degree 3 (monic (x+1000000007)^3)
    assert(polynomial_degree(g) == 3);
    // verify: evaluate gcd at x = -1000000007, should be 0
    double val = poly_evaluate(g, -1000000007.0);
    // this may lose precision in double eval, so check structurally instead:
    // the gcd should be monic and divide both p and q
    Polynomial* rem_p = poly_mod(p, g);
    Polynomial* rem_q = poly_mod(q, g);
    assert(poly_is_zero(rem_p));
    assert(poly_is_zero(rem_q));
    (void)val;

    free_polynomial(rem_p); free_polynomial(rem_q);
    free_polynomial(p); free_polynomial(q); free_polynomial(g);
    free_polynomial(factor1); free_polynomial(factor2);
    rat_clear(&one); rat_clear(&big); rat_clear(&thirteen);
    printf("  test_gmp_poly_gcd_large PASSED\n");
}

void test_gmp_det_large_matrix() {
    // 6x6 matrix with entries up to 100 — det computation involves
    // intermediate products that exceed long long
    // Vandermonde matrix V[i][j] = i^j for i=1..6, j=0..5
    // det = product of (i-j) for all i > j = 1*2*1*3*2*1*4*3*2*1*5*4*3*2*1
    //     = 1! * 2! * 3! * 4! * 5! = 1*2*6*24*120 = 34560
    // (but transposed Vandermonde has same det)
    Matrix* V = matrix_create(6, 6);
    for (int r = 0; r < 6; r++) {
        Rational base = rat_from_int(r + 1);
        Rational power = rat_one_val();
        for (int c = 0; c < 6; c++) {
            rat_set(&MAT_AT(V, r, c), &power);
            Rational tmp; rat_init(&tmp);
            rat_mul(&tmp, &power, &base);
            rat_set(&power, &tmp);
            rat_clear(&tmp);
        }
        rat_clear(&base); rat_clear(&power);
    }

    Rational d; rat_init(&d);
    matrix_det(&d, V);
    // Vandermonde det = prod_{0<=i<j<=5} (x_j - x_i) where x_k = k+1
    // = (2-1)(3-1)(3-2)(4-1)(4-2)(4-3)(5-1)(5-2)(5-3)(5-4)(6-1)(6-2)(6-3)(6-4)(6-5)
    // = 1*2*1*3*2*1*4*3*2*1*5*4*3*2*1 = 34560
    Rational expected = rat_from_int(34560);
    assert(rat_eq(&d, &expected));

    rat_clear(&d); rat_clear(&expected); matrix_free(V);
    printf("  test_gmp_det_large_matrix PASSED\n");
}

void test_gmp_chained_rational_ops() {
    // compute sum 1/1 + 1/2 + 1/3 + ... + 1/100 exactly
    // then verify denominator is lcm(1..100) and the sum is correct
    // by multiplying by lcm and checking we get an integer
    Rational sum = rat_zero_val();
    for (int i = 1; i <= 100; i++) {
        Rational term = rat_from_si(1, i);
        Rational tmp; rat_init(&tmp);
        rat_add(&tmp, &sum, &term);
        rat_set(&sum, &tmp);
        rat_clear(&tmp); rat_clear(&term);
    }
    // sum * 100! should be an integer if computed exactly
    // simpler check: sum is not zero and denominator is huge
    assert(!rat_is_zero(&sum));
    char* s = rat_to_string(&sum);
    // H_100 = 14466636279520351160221518043104131447711/278881500...
    // the string should contain '/' (it's not an integer)
    assert(strstr(s, "/") != NULL);
    // numerator should be very large
    assert(strlen(s) > 20);
    printf("  H_100 = %s\n", s);

    // verify: subtract each 1/k back out, should get zero
    for (int i = 1; i <= 100; i++) {
        Rational term = rat_from_si(1, i);
        Rational tmp; rat_init(&tmp);
        rat_sub(&tmp, &sum, &term);
        rat_set(&sum, &tmp);
        rat_clear(&tmp); rat_clear(&term);
    }
    assert(rat_is_zero(&sum));

    free(s);
    rat_clear(&sum);
    printf("  test_gmp_chained_rational_ops PASSED\n");
}

void test_gmp_factor_large_coefficients() {
    // factor 1000000007x^2 - 1000000007 = 1000000007(x-1)(x+1)
    Polynomial* p = create_polynomial();
    Rational big = rat_from_int(1000000007);
    Rational neg_big; rat_init(&neg_big);
    rat_neg(&neg_big, &big);
    add_term(p, &big, 2);
    add_term(p, &neg_big, 0);
    rat_clear(&big); rat_clear(&neg_big);

    Factorization* f = factorize(p);
    assert(f);
    // content should be 1000000007
    Rational expected_content = rat_from_int(1000000007);
    assert(rat_eq(&f->content, &expected_content));
    rat_clear(&expected_content);
    // should have 2 linear factors
    assert(f->count == 2);
    for (int i = 0; i < 2; i++) {
        assert(polynomial_degree(f->factors[i]) == 1);
        assert(f->multiplicities[i] == 1);
    }

    free_factorization(f); free_polynomial(p);
    printf("  test_gmp_factor_large_coefficients PASSED\n");
}

void test_gmp_inverse_rational_entries() {
    // inverse of [[1/2, 1/3], [1/4, 1/5]]
    // det = 1/10 - 1/12 = 1/60
    // inv = (1/det) * [[1/5, -1/3], [-1/4, 1/2]]
    //     = 60 * [[1/5, -1/3], [-1/4, 1/2]]
    //     = [[12, -20], [-15, 30]]
    Matrix* m = matrix_create(2, 2);
    Rational v;
    v = rat_from_si(1, 2); rat_set(&MAT_AT(m, 0, 0), &v); rat_clear(&v);
    v = rat_from_si(1, 3); rat_set(&MAT_AT(m, 0, 1), &v); rat_clear(&v);
    v = rat_from_si(1, 4); rat_set(&MAT_AT(m, 1, 0), &v); rat_clear(&v);
    v = rat_from_si(1, 5); rat_set(&MAT_AT(m, 1, 1), &v); rat_clear(&v);

    Matrix* inv = matrix_inverse(m);
    assert(inv);

    Rational e;
    e = R(12); assert(rat_eq(&MAT_AT(inv, 0, 0), &e)); rat_clear(&e);
    e = R(-20); assert(rat_eq(&MAT_AT(inv, 0, 1), &e)); rat_clear(&e);
    e = R(-15); assert(rat_eq(&MAT_AT(inv, 1, 0), &e)); rat_clear(&e);
    e = R(30); assert(rat_eq(&MAT_AT(inv, 1, 1), &e)); rat_clear(&e);

    matrix_free(m); matrix_free(inv);
    printf("  test_gmp_inverse_rational_entries PASSED\n");
}

/***********************************************************************
*  run all tests
*************************************************************************/

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

    printf("\nComplex / Stress Tests:\n");
    test_parse_expand_product();
    test_expand_then_factor_roundtrip();
    test_differentiate_integrate_inverse();
    test_high_degree_polynomial_multiply();
    test_solve_quartic();
    test_solve_degree5_partial();
    test_parse_complex_expression();
    test_rational_arithmetic_stress();
    test_power_of_binomial();
    test_gmp_large_coefficients();

    printf("\nFormatter Tests:\n");
    test_fmt_add_negative_rhs();
    test_fmt_mul_coeff_one();
    test_fmt_mul_coeff_neg_one();
    test_fmt_neg_coeff_in_sum();
    test_fmt_pow_neg_one();
    test_fmt_pow_neg_two();
    test_fmt_div_parens();
    test_fmt_pow_poly_base();
    test_fmt_chained_add_negatives();
    test_fmt_integration_output();

    printf("\nGMP Precision Tests:\n");
    test_gmp_factorial_fraction();
    test_gmp_hilbert_5x5();
    test_gmp_poly_gcd_large();
    test_gmp_det_large_matrix();
    test_gmp_chained_rational_ops();
    test_gmp_factor_large_coefficients();
    test_gmp_inverse_rational_entries();

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
