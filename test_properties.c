// Regression tests for every defect found by probing the CLI, plus seeded
// property tests that check invariants on generated inputs instead of fixed outputs:
//   factor:     content · Π f_i^m_i == input, and known irreducibles come back exactly
//   roots:      every isolating interval brackets a sign change; real-root counts match
//   integrate:  d/dx F == f exactly, on integrands built backwards from known pieces
//   linear:     A · solve(A, A·x) == A·x exactly, det(A)·det(A^-1) == 1
//   printing:   parse(print(p)) == p for large random polynomials
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include "expression.h"
#include "factor.h"
#include "solver.h"
#include "matrix.h"

#define SEED 20261113ULL
#define TRIALS 60

// --- deterministic PRNG (xorshift64*) ---

static unsigned long long rng_state = SEED;

static long rnd(long lo, long hi) { // inclusive range
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    unsigned long long v = rng_state * 2685821657736338717ULL;
    return lo + (long)(v % (unsigned long long)(hi - lo + 1));
}

// --- helpers ---

// run a CLI-style action on an input string and return the printed result (caller frees)
static char* simplify_str(const char* in) {
    Expr* e = parse_and_simplify(in);
    assert(e);
    char* s = expr_to_string(e);
    free_expr(e);
    return s;
}

static char* factor_str(const char* in) {
    Expr* e = parse_and_simplify(in);
    assert(e && e->type == EXPR_POLY);
    Factorization* f = factorize(e->poly);
    char* s = factorization_to_string(f);
    free_factorization(f);
    free_expr(e);
    return s;
}

static Expr* integrate_str(const char* in) {
    Expr* e = parse_and_simplify(in);
    assert(e);
    Expr* r = integrate(e, "x");
    free_expr(e);
    return r;
}

static void expect_str(char* got, const char* want) {
    if (strcmp(got, want) != 0) {
        fprintf(stderr, "    expected: %s\n    got:      %s\n", want, got);
        assert(0);
    }
    free(got);
}

static Polynomial* poly_si(const long* c, int deg) { // c[i] is the coefficient of x^i
    Polynomial* p = create_polynomial();
    for (int i = deg; i >= 0; i--) {
        Rational r = rat_from_int(c[i]);
        add_term(p, &r, i);
        rat_clear(&r);
    }
    return p;
}

static bool poly_equal(const Polynomial* a, const Polynomial* b) {
    Polynomial* d = poly_subtract(a, b);
    bool eq = poly_is_zero(d);
    free_polynomial(d);
    return eq;
}

/***********************************************************************
*  regression tests: one per defect found by running the CLI
*************************************************************************/

static void test_reg_big_literals_exact(void) {
    expect_str(simplify_str("x^2 - 10000000000000000000000*x"), "x^2 - 10000000000000000000000x");
    expect_str(simplify_str("9007199254740993x"), "9007199254740993x"); // 2^53 + 1
    expect_str(simplify_str("0.1x + 0.2x"), "3/10x");
    expect_str(simplify_str("1.5e3x + 2e-2"), "1500x + 1/50");
    printf("  test_reg_big_literals_exact PASSED\n");
}

static void test_reg_matrix_big_entries(void) {
    Matrix* m = matrix_parse("[[10000000000000000000000,0],[0,1]]");
    assert(m);
    Rational d;
    rat_init(&d);
    matrix_det(&d, m);
    expect_str(rat_to_string(&d), "10000000000000000000000");
    rat_clear(&d);
    matrix_free(m);
    printf("  test_reg_matrix_big_entries PASSED\n");
}

static void test_reg_output_not_truncated(void) {
    Expr* e = parse_and_simplify("(x+1)^500");
    assert(e && e->type == EXPR_POLY);
    char* s = expr_to_string(e);
    mpz_t binom;
    mpz_init(binom);
    mpz_bin_uiui(binom, 500, 250); // 150 digits, past the old 128-byte buffer
    char* want = mpz_get_str(NULL, 10, binom);
    size_t cap = strlen(want) + 8;
    char* needle = malloc(cap);
    snprintf(needle, cap, " %sx^250 ", want);
    assert(strstr(s, needle) != NULL);
    free(needle); free(want); free(s);
    mpz_clear(binom);
    free_expr(e);
    printf("  test_reg_output_not_truncated PASSED\n");
}

static void test_reg_rejects_non_polynomials(void) {
    assert(parse_and_simplify("1/0") == NULL);
    assert(parse_and_simplify("x^(1/2)") == NULL);
    assert(parse_and_simplify("x^1000000000") == NULL);
    printf("  test_reg_rejects_non_polynomials PASSED\n");
}

// trailing or malformed input is an error, never a silently shortened expression
static void test_reg_parser_rejects_partial(void) {
    const char* bad[] = { "x^2 = 4", "2x $ 3", "2x)", "(2x", "2 +", "*3" };
    for (size_t i = 0; i < sizeof bad / sizeof bad[0]; i++) assert(parse_and_simplify(bad[i]) == NULL);
    expect_str(simplify_str("  3x  "), "3x");
    expect_str(simplify_str("-x^2 + 1"), "-x^2 + 1");  // used to parse as (-x)^2
    expect_str(simplify_str("-2^2"), "-4");
    expect_str(simplify_str("x^-1 * x^2"), "x");
    expect_str(simplify_str("-2*-3x"), "6x");
    printf("  test_reg_parser_rejects_partial PASSED\n");
}

static void test_reg_factor_content_and_big_roots(void) {
    expect_str(factor_str("2x^2 - 3x + 1"), "(2x - 1)(x - 1)");       // the 2 used to be dropped
    expect_str(factor_str("6x^2+5x+1"), "(2x + 1)(3x + 1)");
    expect_str(factor_str("(x - 4294967296)^3"), "(x - 4294967296)^3"); // constant -2^96
    expect_str(factor_str("x^2 - 10000000000000000000000*x"), "(x)(x - 10000000000000000000000)");
    printf("  test_reg_factor_content_and_big_roots PASSED\n");
}

static void test_reg_factor_beyond_rational_roots(void) {
    expect_str(factor_str("x^4 - 5x^2 + 6"), "(x^2 - 3)(x^2 - 2)");
    expect_str(factor_str("x^4+4"), "(x^2 - 2x + 2)(x^2 + 2x + 2)");
    expect_str(factor_str("x^4+1"), "(x^4 + 1)");                 // splits mod every prime
    expect_str(factor_str("x^4-10x^2+1"), "(x^4 - 10x^2 + 1)");   // splits mod every prime
    expect_str(factor_str("(x^2-2)^2"), "(x^2 - 2)^2");           // square-free decomposition
    expect_str(factor_str("(2x^2-1)(3x^2+1)"), "(2x^2 - 1)(3x^2 + 1)");
    expect_str(factor_str("x^8 + 98x^4 + 1"), "(x^4 - 4x^3 + 8x^2 + 4x + 1)(x^4 + 4x^3 + 8x^2 - 4x + 1)");
    printf("  test_reg_factor_beyond_rational_roots PASSED\n");
}

// x^n - 1 has one cyclotomic factor per divisor of n
static void test_reg_cyclotomic_counts(void) {
    int ns[] = { 6, 12, 30, 64 }, divisors[] = { 4, 6, 8, 7 };
    for (int i = 0; i < 4; i++) {
        char in[32];
        snprintf(in, sizeof in, "x^%d - 1", ns[i]);
        Expr* e = parse_and_simplify(in);
        Factorization* f = factorize(e->poly);
        assert(f->count == divisors[i]);
        free_factorization(f);
        free_expr(e);
    }
    printf("  test_reg_cyclotomic_counts PASSED\n");
}

static void test_reg_solve_irrational_exact(void) {
    Expr* e = parse_and_simplify("x^3 - 2");
    Solutions* s = solve_polynomial(e->poly);
    assert(s->num_rational == 0 && s->num_irr == 1 && s->irr[0].num_real == 1);
    double lo = rat_to_double(&s->irr[0].real[0].lo), hi = rat_to_double(&s->irr[0].real[0].hi);
    assert(lo <= cbrt(2.0) && cbrt(2.0) <= hi && hi - lo < 1e-11);
    char* str = solutions_to_string(s);
    assert(strstr(str, "x ~ 1.2599210498") && strstr(str, "2 non-real roots"));
    free(str);
    free_solutions(s);
    free_expr(e);
    Expr* q = parse_and_simplify("x^2 + x + 1");
    s = solve_polynomial(q->poly);
    expect_str(solutions_to_string(s), "x = -1/2 +/- 1/2*sqrt(3)*i");
    free_solutions(s);
    free_expr(q);
    printf("  test_reg_solve_irrational_exact PASSED\n");
}

static void test_reg_integrate_complete(void) {
    Expr* r = integrate_str("1/((x-1)(x^2+1))"); // used to silently drop two terms
    expect_str(expr_to_string(r), "1/2 * ln|x - 1| - 1/4 * ln|x^2 + 1| - 1/2 * atan(x)");
    free_expr(r);
    r = integrate_str("1/(2x^2 - 3x + 1)");      // used to come out doubled
    expect_str(expr_to_string(r), "-ln|2x - 1| + ln|x - 1|");
    free_expr(r);
    r = integrate_str("1/(x^2+1)^2");
    expect_str(expr_to_string(r), "1/2x / (x^2 + 1) + 1/2 * atan(x)");
    free_expr(r);
    assert(integrate_str("1/(x^3-2)") == NULL && integrate_error() != NULL); // fails loudly
    printf("  test_reg_integrate_complete PASSED\n");
}

/***********************************************************************
*  property tests on generated inputs
*************************************************************************/

static long gcd_l(long a, long b) { return b == 0 ? a : gcd_l(b, a % b); }

// a random factor whose irreducibility over Q is known by construction:
// linear a·x + b, x^2 + k (k > 0), or x^2 - d (d not a square); returns real-root count
static Polynomial* known_irreducible(int* real_roots) {
    long c[3] = { 0, 0, 0 };
    int kind = (int)rnd(0, 2);
    if (kind == 0) { // primitive a·x + b, i.e. gcd(a, |b|) == 1
        do { c[1] = rnd(1, 9); c[0] = rnd(-12, 12); } while (gcd_l(c[1], labs(c[0])) != 1);
        *real_roots = 1;
        return poly_si(c, 1);
    }
    c[2] = 1;
    if (kind == 1) { c[0] = rnd(1, 30); *real_roots = 0; }
    else {
        static const long nonsquares[] = { 2, 3, 5, 6, 7, 10, 11, 13, 15, 17 };
        c[0] = -nonsquares[rnd(0, 9)];
        *real_roots = 2;
    }
    return poly_si(c, 2);
}

// true if p was already generated (keeps the factor set distinct)
static bool already_used(Polynomial** fs, int n, const Polynomial* p) {
    for (int i = 0; i < n; i++) if (poly_equal(fs[i], p)) return true;
    return false;
}

typedef struct {
    Polynomial* product;
    Polynomial* factors[4];
    int mult[4];
    int real_roots[4];
    int count;
} Generated;

static Generated generate_product(void) {
    Generated g = { NULL, { NULL }, { 0 }, { 0 }, (int)rnd(1, 4) };
    Rational content = rat_from_si(rnd(1, 7) * (rnd(0, 1) ? 1 : -1), rnd(1, 5));
    g.product = poly_constant(&content);
    for (int i = 0; i < g.count; i++) {
        Polynomial* f;
        do { f = known_irreducible(&g.real_roots[i]); }
        while (already_used(g.factors, i, f) && (free_polynomial(f), 1));
        g.factors[i] = f;
        g.mult[i] = (int)rnd(1, 3);
        Polynomial* pw = poly_pow(f, (unsigned)g.mult[i]);
        Polynomial* t = poly_multiply(g.product, pw);
        free_polynomial(g.product); free_polynomial(pw);
        g.product = t;
    }
    rat_clear(&content);
    return g;
}

static void free_generated(Generated* g) {
    free_polynomial(g->product);
    for (int i = 0; i < g->count; i++) free_polynomial(g->factors[i]);
}

static void test_prop_factor_roundtrip(void) {
    for (int trial = 0; trial < TRIALS; trial++) {
        Generated g = generate_product();
        Factorization* f = factorize(g.product);
        Polynomial* back = poly_constant(&f->content);
        for (int i = 0; i < f->count; i++) {
            Polynomial* pw = poly_pow(f->factors[i], (unsigned)f->multiplicities[i]);
            Polynomial* t = poly_multiply(back, pw);
            free_polynomial(back); free_polynomial(pw);
            back = t;
        }
        assert(poly_equal(back, g.product));            // exact reconstruction
        assert(f->count == g.count);                     // known irreducibles: same multiset
        for (int i = 0; i < g.count; i++) {
            bool found = false;
            for (int j = 0; j < f->count; j++)
                if (poly_equal(f->factors[j], g.factors[i]) && f->multiplicities[j] == g.mult[i]) found = true;
            assert(found);
        }
        free_polynomial(back);
        free_factorization(f);
        free_generated(&g);
    }
    printf("  test_prop_factor_roundtrip PASSED (%d trials)\n", TRIALS);
}

static void test_prop_root_isolation(void) {
    for (int trial = 0; trial < TRIALS; trial++) {
        Generated g = generate_product();
        Solutions* s = solve_polynomial(g.product);
        int expect_real = 0, got_real = s->num_rational;
        for (int i = 0; i < g.count; i++) expect_real += g.real_roots[i] * g.mult[i];
        Rational v;
        rat_init(&v);
        for (int i = 0; i < s->num_rational; i++) { // rational roots are exact zeros
            poly_eval_rat(&v, g.product, &s->rational_roots[i]);
            assert(rat_is_zero(&v));
        }
        for (int i = 0; i < s->num_irr; i++) {
            const IrrationalFactor* e = &s->irr[i];
            got_real += e->num_real * e->multiplicity;
            for (int k = 0; k < e->num_real; k++) { // each interval brackets a sign change
                Rational a;
                rat_init(&a);
                poly_eval_rat(&a, e->factor, &e->real[k].lo);
                poly_eval_rat(&v, e->factor, &e->real[k].hi);
                assert(mpq_sgn(a.val) * mpq_sgn(v.val) < 0);
                rat_clear(&a);
            }
        }
        assert(got_real == expect_real);
        rat_clear(&v);
        free_solutions(s);
        free_generated(&g);
    }
    printf("  test_prop_root_isolation PASSED (%d trials)\n", TRIALS);
}

// --- exact derivative of an antiderivative produced by integrate() ---
// Every term shape the integrator emits has a rational derivative over Q:
//   c·ln|P|                        -> c·P'/P
//   c·atan(L)                      -> c·L'/(1 + L^2)
//   c·sqrt(m)·atan(L/sqrt(m))      -> c·m·L'/(L^2 + m)
//   c·sqrt(m)·ln|(L-√m)/(L+√m)|    -> c·2m·L'/(L^2 - m)
//   N/D, -1·(N/D), polynomial      -> quotient rule
typedef struct { Polynomial* num; Polynomial* den; } RF;

static RF rf_make(Polynomial* n, Polynomial* d) { return (RF){ n, d }; }

static void rf_free(RF* r) { free_polynomial(r->num); free_polynomial(r->den); }

static void rf_add_into(RF* acc, RF t) { // consumes t
    Polynomial *a = poly_multiply(acc->num, t.den), *b = poly_multiply(t.num, acc->den);
    Polynomial *n = poly_add(a, b), *d = poly_multiply(acc->den, t.den);
    free_polynomial(a); free_polynomial(b); rf_free(acc); rf_free(&t);
    *acc = rf_make(n, d);
}

static void rf_scale(RF* r, const Rational* c) {
    Polynomial* n = poly_scale(r->num, c);
    free_polynomial(r->num);
    r->num = n;
}

static Polynomial* poly_one_t(void) { Rational o = rat_one_val(); Polynomial* p = poly_constant(&o); rat_clear(&o); return p; }

// (N/D)' = (N'D - ND')/D^2
static RF rf_quotient_rule(const Polynomial* N, const Polynomial* D) {
    Polynomial *np = poly_derivative(N), *dp = poly_derivative(D);
    Polynomial *a = poly_multiply(np, D), *b = poly_multiply(N, dp);
    RF r = rf_make(poly_subtract(a, b), poly_multiply(D, D));
    free_polynomial(np); free_polynomial(dp); free_polynomial(a); free_polynomial(b);
    return r;
}

// k·L'/(L^2 + s) with s = shift (m, -m, or 1)
static RF rf_arc(const Polynomial* L, const Rational* k, const Rational* shift) {
    Polynomial *lp = poly_derivative(L), *l2 = poly_multiply(L, L), *sh = poly_constant(shift);
    RF r = rf_make(poly_scale(lp, k), poly_add(l2, sh));
    free_polynomial(lp); free_polynomial(l2); free_polynomial(sh);
    return r;
}

static const Rational* sqrt_arg(const Expr* e) { // the m in sqrt(m), or NULL
    return (e->type == EXPR_SQRT && e->fn_arg->type == EXPR_CONST) ? &e->fn_arg->constant : NULL;
}

// derivative of sqrt(m)·atan(L/sqrt m) or sqrt(m)·ln|(L-√m)/(L+√m)|; NULL shape -> false
static bool rf_sqrt_term(const Expr* rest, const Rational* m, RF* out) {
    Rational k, sh;
    rat_init(&k); rat_init(&sh);
    bool ok = false;
    if (rest->type == EXPR_ATAN && rest->fn_arg->type == EXPR_BINOP && rest->fn_arg->binop.op == OP_DIV) {
        rat_set(&k, m); rat_set(&sh, m);                         // m·L'/(L^2 + m)
        *out = rf_arc(rest->fn_arg->binop.left->poly, &k, &sh);
        ok = true;
    } else if (rest->type == EXPR_LN && rest->ln_arg->type == EXPR_BINOP) {
        rat_add(&k, m, m); rat_neg(&sh, m);                      // 2m·L'/(L^2 - m)
        *out = rf_arc(rest->ln_arg->binop.left->binop.left->poly, &k, &sh);
        ok = true;
    }
    rat_clear(&k); rat_clear(&sh);
    return ok;
}

// derivative of one emitted term c·(...)
static bool rf_scaled_term(const Expr* t, RF* out) {
    const Rational* c = &t->binop.left->constant;
    const Expr* g = t->binop.right;
    bool ok = true;
    if (g->type == EXPR_LN && g->ln_arg->type == EXPR_POLY) {
        *out = rf_make(poly_derivative(g->ln_arg->poly), poly_copy(g->ln_arg->poly));
    } else if (g->type == EXPR_ATAN && g->fn_arg->type == EXPR_POLY) {
        Rational one = rat_one_val();
        *out = rf_arc(g->fn_arg->poly, &one, &one);
        rat_clear(&one);
    } else if (g->type == EXPR_BINOP && g->binop.op == OP_DIV) {
        *out = rf_quotient_rule(g->binop.left->poly, g->binop.right->poly);
    } else if (g->type == EXPR_BINOP && g->binop.op == OP_MUL && sqrt_arg(g->binop.left)) {
        ok = rf_sqrt_term(g->binop.right, sqrt_arg(g->binop.left), out);
    } else {
        ok = false;
    }
    if (ok) rf_scale(out, c);
    return ok;
}

// exact derivative of F, accumulated into acc; false on an unrecognized shape
static bool rf_derivative(const Expr* F, RF* acc) {
    RF t;
    if (F->type == EXPR_BINOP && F->binop.op == OP_ADD)
        return rf_derivative(F->binop.left, acc) && rf_derivative(F->binop.right, acc);
    if (F->type == EXPR_POLY) t = rf_make(poly_derivative(F->poly), poly_one_t());
    else if (F->type == EXPR_BINOP && F->binop.op == OP_DIV) t = rf_quotient_rule(F->binop.left->poly, F->binop.right->poly);
    else if (F->type == EXPR_BINOP && F->binop.op == OP_MUL && F->binop.left->type == EXPR_CONST) {
        if (!rf_scaled_term(F, &t)) return false;
    } else return false;
    rf_add_into(acc, t);
    return true;
}

// F' == num/den exactly
static bool derivative_is(const Expr* F, const Polynomial* num, const Polynomial* den) {
    RF acc = rf_make(create_polynomial(), poly_one_t());
    bool ok = rf_derivative(F, &acc);
    if (ok) {
        Polynomial *x = poly_multiply(acc.num, den), *y = poly_multiply(num, acc.den);
        ok = poly_equal(x, y);
        free_polynomial(x); free_polynomial(y);
    }
    rf_free(&acc);
    return ok;
}

// backward generation: integrands whose antiderivatives are known to be elementary over Q
static void test_prop_integrate_derivative(void) {
    for (int trial = 0; trial < TRIALS; trial++) {
        Generated g = generate_product();
        long nc[6];
        int nd = (int)rnd(0, polynomial_degree(g.product) + 1);
        if (nd > 5) nd = 5;
        for (int i = 0; i <= nd; i++) nc[i] = rnd(-9, 9);
        nc[nd] = nc[nd] ? nc[nd] : 1;
        Polynomial* num = poly_si(nc, nd);
        Expr* integrand = create_div(create_poly_expr(poly_copy(num)), create_poly_expr(poly_copy(g.product)));
        Expr* F = integrate(integrand, "x");
        if (!F) fprintf(stderr, "    integrate failed: %s\n", integrate_error());
        assert(F != NULL);
        if (!derivative_is(F, num, g.product)) {
            char *fs = expr_to_string(integrand), *Fs = expr_to_string(F);
            fprintf(stderr, "    d/dx F != f exactly\n    f = %s\n    F = %s\n", fs, Fs);
            free(fs); free(Fs);
            assert(0);
        }
        free_expr(F);
        free_expr(integrand);
        free_polynomial(num);
        free_generated(&g);
    }
    printf("  test_prop_integrate_derivative PASSED (%d trials)\n", TRIALS);
}

static Matrix* random_matrix(int n) {
    Matrix* A = matrix_create(n, n);
    for (int i = 0; i < n * n; i++) rat_set_si(&A->data[i], rnd(-20, 20), rnd(1, 4));
    return A;
}

static void test_prop_linear_systems(void) {
    for (int trial = 0; trial < TRIALS; trial++) {
        int n = (int)rnd(1, 6);
        Matrix* A = random_matrix(n);
        Rational det;
        rat_init(&det);
        matrix_det(&det, A);
        if (rat_is_zero(&det)) { rat_clear(&det); matrix_free(A); continue; }
        Vector* x = vector_create(n);
        Vector* b = vector_create(n);
        Rational t;
        rat_init(&t);
        for (int i = 0; i < n; i++) rat_set_si(&x->data[i], rnd(-50, 50), rnd(1, 7));
        for (int r = 0; r < n; r++) // b = A·x
            for (int c = 0; c < n; c++) {
                rat_mul(&t, &MAT_AT(A, r, c), &x->data[c]);
                rat_add(&b->data[r], &b->data[r], &t);
            }
        LinSysResult* res = linsys_solve(A, b);
        assert(res->status == LINSYS_UNIQUE);
        for (int i = 0; i < n; i++) assert(rat_eq(&res->solution->data[i], &x->data[i]));
        Matrix* inv = matrix_inverse(A);
        matrix_det(&t, inv);
        rat_mul(&t, &t, &det);
        assert(rat_is_one(&t)); // det(A)·det(A^-1) == 1
        matrix_free(inv);
        linsys_result_free(res);
        rat_clear(&t); rat_clear(&det);
        vector_free(x); vector_free(b); matrix_free(A);
    }
    printf("  test_prop_linear_systems PASSED (%d trials)\n", TRIALS);
}

// parse(print(p)) == p, with coefficients far beyond 64 bits
static void test_prop_print_parse_roundtrip(void) {
    for (int trial = 0; trial < TRIALS; trial++) {
        Polynomial* p = create_polynomial();
        for (int e = (int)rnd(0, 12); e >= 0; e--) {
            Rational c;
            rat_init(&c);
            mpz_ui_pow_ui(mpq_numref(c.val), (unsigned long)rnd(2, 97), (unsigned long)rnd(1, 60));
            mpz_set_ui(mpq_denref(c.val), (unsigned long)rnd(1, 30));
            mpq_canonicalize(c.val);
            if (rnd(0, 1)) rat_neg(&c, &c);
            add_term(p, &c, e);
            rat_clear(&c);
        }
        char* s = poly_to_string(p);
        Expr* back = parse_and_simplify(s);
        assert(back && back->type == EXPR_POLY && poly_equal(back->poly, p));
        free_expr(back);
        free(s);
        free_polynomial(p);
    }
    printf("  test_prop_print_parse_roundtrip PASSED (%d trials)\n", TRIALS);
}

void run_property_tests(void) {
    const char* env = getenv("POLY_SEED"); // optional override to sweep other seeds
    if (env) rng_state = strtoull(env, NULL, 10);
    printf("\nRegression Tests (CLI probes):\n");
    test_reg_big_literals_exact();
    test_reg_matrix_big_entries();
    test_reg_output_not_truncated();
    test_reg_rejects_non_polynomials();
    test_reg_parser_rejects_partial();
    test_reg_factor_content_and_big_roots();
    test_reg_factor_beyond_rational_roots();
    test_reg_cyclotomic_counts();
    test_reg_solve_irrational_exact();
    test_reg_integrate_complete();
    printf("\nProperty Tests (seed %s):\n", env ? env : "default");
    test_prop_factor_roundtrip();
    test_prop_root_isolation();
    test_prop_integrate_derivative();
    test_prop_linear_systems();
    test_prop_print_parse_roundtrip();
}
