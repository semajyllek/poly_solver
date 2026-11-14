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
#include "logpart.h"
#include "test_subres_ref.h"
#include "rootform.h"
#include "definite.h"
#include "constfold.h"
#include "qsurd.h"
#include "modgcd.h"
#include "risch.h"
#include "expfield.h"
#include <complex.h>

static bool derivative_is(const Expr* F, const Polynomial* num, const Polynomial* den);

typedef struct { double x0, h; } Stencil;

// 5-point central difference of F at st.x0 with step st.h
static double stencil(const Expr* F, Stencil st, bool* ok) {
    double x0 = st.x0, h = st.h;
    return (-expr_eval(F, x0 + 2 * h, ok) + 8 * expr_eval(F, x0 + h, ok)
            - 8 * expr_eval(F, x0 - h, ok) + expr_eval(F, x0 - 2 * h, ok)) / (12 * h);
}

// F' == num/den in double precision at sample points where the difference quotient has
// converged (two step sizes agree, i.e. no pole nearby); for answers with radicals and trig
// constants, which the exact recognizer does not model
static bool numeric_derivative_is(const Expr* F, const Polynomial* num, const Polynomial* den) {
    const double xs[] = { -2.7, -1.3, -0.45, 0.35, 1.15, 2.45, -3.9, 3.7 };
    int checked = 0;
    for (size_t i = 0; i < sizeof xs / sizeof xs[0]; i++) {
        double x0 = xs[i], d = poly_evaluate(den, x0);
        bool ok = d != 0;
        double f = poly_evaluate(num, x0) / d, d1 = stencil(F, (Stencil){ x0, 1e-3 }, &ok),
               d2 = stencil(F, (Stencil){ x0, 5e-4 }, &ok);
        if (!ok || fabs(d1 - d2) > 1e-9 * (1 + fabs(f))) continue; // not converged: pole nearby
        if (fabs(d2 - f) > 1e-7 * (1 + fabs(f))) return false;
        checked++;
    }
    return checked >= 3;
}

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

typedef Expr* (*IntegrateFn)(Expr*, const char*);

static Expr* integrate_with(IntegrateFn fn, const char* in) {
    Expr* e = parse_and_simplify(in);
    assert(e);
    Expr* r = fn(e, "x");
    free_expr(e);
    return r;
}

static Expr* integrate_str(const char* in) { return integrate_with(integrate, in); }

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
// Lazard–Rioboo–Trager: irreducible denominator factors of degree >= 3 (RootSum kept)
static void test_reg_integrate_lrt(void) {
    struct { const char* f; const char* F; } cases[] = {
        { "1/(x^3-2)", "RootSum(108t^3 - 1, t -> t * ln(x - 6t))" },        // used to be refused
        { "1/(x^4+1)", "RootSum(256t^4 + 1, t -> t * ln(x + 4t))" },
        { "x/(x^4+1)", "1/2 * atan(x^2)" },                                 // complex pair, LogToAtan
        { "x^2/(x^6+1)", "1/3 * atan(x^3)" },                               // like terms merged
        { "(x^2+1)/(x^4-x^2+1)", "atan(x^3) + atan(x)" },                   // continuous atans
        { "x/(x^4-2)", "-1/8 * sqrt(2) * ln(abs((x^2 + sqrt(2)) / (x^2 - sqrt(2))))" }, // real pair
        { "3x^2/(x^3-2)", "ln(abs(x^3 - 2))" },                                  // rational log coefficient
        { "x/(x^4+3)", "1/6 * sqrt(3) * atan(x^2 / sqrt(3))" },             // complex pair, irrational i-part
        { "x^3/(x^8+5)", "1/20 * sqrt(5) * atan(x^4 / sqrt(5))" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        Expr* r = integrate_with(integrate_rootsum, cases[i].f);
        assert(r != NULL);
        expect_str(expr_to_string(r), cases[i].F);
        free_expr(r);
    }
    printf("  test_reg_integrate_lrt PASSED\n");
}

// RootSum expanded into radicals, cos/sin(kπ/N) and real ln/atan when the roots are known
// (binomials and rescaled cyclotomics); every string was checked to 60 digits with SymPy
static void test_reg_integrate_radicals(void) {
    struct { const char* f; const char* F; } cases[] = {
        { "1/(x^3-2)", // real root + one pair
          "1/6 * 2^(1/3) * ln(abs(x - 2^(1/3))) - 1/12 * 2^(1/3) * ln(abs(x^2 + 2^(1/3) * x + 2^(2/3))) - 1/6 * 2^(1/3) * sqrt(3) * atan(1/3 * 2^(2/3) * sqrt(3) * x + 1/3 * sqrt(3))" },
        { "1/(x^4+1)", // textbook sqrt(2) form
          "-1/8 * sqrt(2) * ln(abs(x^2 - sqrt(2) * x + 1)) + 1/4 * sqrt(2) * atan(sqrt(2) * x - 1) + 1/8 * sqrt(2) * ln(abs(x^2 + sqrt(2) * x + 1)) + 1/4 * sqrt(2) * atan(sqrt(2) * x + 1)" },
        { "1/(x^6+1)", // cyclotomic factors, table angles
          "1/3 * atan(x) - 1/12 * sqrt(3) * ln(abs(x^2 - sqrt(3) * x + 1)) + 1/6 * atan(2 * x - sqrt(3)) + 1/12 * sqrt(3) * ln(abs(x^2 + sqrt(3) * x + 1)) + 1/6 * atan(2 * x + sqrt(3))" },
        { "1/(x^5-1)", // cyclotomic, symbolic cos(2*pi/5)
          "1/5 * ln(abs(x - 1)) + 1/5 * cos(2*pi/5) * ln(abs(x^2 - 2 * cos(2*pi/5) * x + 1)) - 2/5 * sin(2*pi/5) * atan(csc(2*pi/5) * x - cot(2*pi/5)) - 1/5 * cos(pi/5) * ln(abs(x^2 + 2 * cos(pi/5) * x + 1)) - 2/5 * sin(pi/5) * atan(csc(pi/5) * x + cot(pi/5))" },
        { "1/(x^6-8)", // Q(x^2) = 4·Φ_3(x^2/2)
          "1/48 * sqrt(2) * ln(abs((x - sqrt(2)) / (x + sqrt(2)))) + 1/96 * sqrt(2) * ln(abs(x^2 - sqrt(2) * x + 2)) - 1/48 * sqrt(6) * atan(1/3 * sqrt(6) * x - 1/3 * sqrt(3)) - 1/96 * sqrt(2) * ln(abs(x^2 + sqrt(2) * x + 2)) - 1/48 * sqrt(6) * atan(1/3 * sqrt(6) * x + 1/3 * sqrt(3))" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        Expr* f = parse_and_simplify(cases[i].f);
        Expr* r = integrate_str(cases[i].f);
        assert(r != NULL && f && f->type == EXPR_RATIONAL_FN);
        assert(numeric_derivative_is(r, f->ratfn.num, f->ratfn.den));
        expect_str(expr_to_string(r), cases[i].F);
        free_expr(r); free_expr(f);
    }
    printf("  test_reg_integrate_radicals PASSED (%zu cases)\n", sizeof cases / sizeof cases[0]);
}

// definite integrals: exact values (folded constants, atan jump corrections, a pole of an
// atan argument at a bound), divergence on a pole in [a, b], and the RootSum fallback
// "inf" / "-inf" → ±1, else 0 with *r set to the rational
static int bound_str(Rational* r, const char* text) {
    if (strcmp(text, "inf") == 0 || strcmp(text, "-inf") == 0) return text[0] == '-' ? -1 : 1;
    Expr* e = parse_and_simplify(text);
    poly_coeff(r, e->poly, 0);
    free_expr(e);
    return 0;
}

static DefiniteResult definite_str(const char* f, const char* a, const char* b) {
    Expr* e = parse_and_simplify(f);
    Rational ra, rb;
    rat_init(&ra); rat_init(&rb);
    Bounds ab = { { &ra, bound_str(&ra, a) }, { &rb, bound_str(&rb, b) } };
    DefiniteResult r = integrate_definite(e, &ab);
    rat_clear(&ra); rat_clear(&rb);
    free_expr(e);
    return r;
}

static void test_reg_integrate_definite(void) {
    const char* cases[][4] = {
        { "1/(x+1)", "0", "1", "ln(2)" },
        { "1/(x^2+1)", "1", "0", "-1/4 * pi" },
        { "1/((x+1)(x+2))", "0", "1", "ln(4/3)" },
        { "1/(x^2+x+1)", "0", "1", "1/9 * sqrt(3) * pi" },
        { "1/(x^4+x^2+1)", "-1", "1", "1/6 * sqrt(3) * pi + 1/2 * ln(3)" },
        { "1/(x^2+1)^2", "0", "1", "1/8 * pi + 1/4" },
        { "(x^2+1)/(x^4+3x^2+1)", "-1", "1", "1/5 * sqrt(5) * pi" },     // atan(N/(sqrt(5)·x)) jumps at 0
        { "(x^2+1)/(x^4+4x^2+1)", "0", "1", "1/12 * sqrt(6) * pi" },     // that pole at a bound
        { "1/(x^2-2)", "-1", "1", "-sqrt(2) * ln(1 + sqrt(2))" },             // ln|1-√2| - ln|1+√2| = -2·ln(1+√2)
        { "1/(x^4+1)", "0", "1", "1/4 * sqrt(2) * ln(1 + sqrt(2)) + 1/8 * sqrt(2) * pi" }, // atan(1±√2) pair → π/2
        { "1/(x^4+1)", "-1", "0", "1/4 * sqrt(2) * ln(1 + sqrt(2)) + 1/8 * sqrt(2) * pi" }, // pair with u < 0: -π/2
        { "1/(x^6+1)", "0", "1", "1/6 * sqrt(3) * ln(2 + sqrt(3)) + 1/6 * pi" },           // 7+4√3 = (2+√3)²
        { "(x+1)/(x^4+1)", "0", "1", "1/4 * sqrt(2) * ln(1 + sqrt(2)) + 1/8 * sqrt(2) * pi + 1/8 * pi" },
        { "1/(x^2-8)", "0", "2", "-1/4 * sqrt(2) * ln(1 + sqrt(2))" },        // |2-2√2| = 2·√2·(√2-1)
        { "1/(x^4+1)", "-2", "-1", "-1/4 * sqrt(2) * ln(1 + 3 * sqrt(2)) + 1/8 * sqrt(2) * ln(17) - 1/4 * sqrt(2) * atan(2/3 * sqrt(2)) + 1/8 * sqrt(2) * pi" }, // ln(-1+3√2) = ln 17 - ln(1+3√2)
        { "1/(x^4+1)", "-3", "-1/2", "-1/4 * sqrt(2) * ln(5 + 19 * sqrt(2)) + 1/8 * sqrt(2) * ln(697) - 1/4 * sqrt(2) * atan(25/12 * sqrt(2)) + 1/4 * sqrt(2) * pi" }, // ln(697²)/16
        { "1/(x^2-3)", "2", "3", "1/3 * sqrt(3) * ln(1 + sqrt(3)) - 1/6 * sqrt(3) * ln(2)" },   // ln 4 = 2 ln 2
        { "x^2", "0", "3", "9" },
        { "x", "2", "2", "0" },
        { "1/(x^4+1)", "-inf", "inf", "1/2 * sqrt(2) * pi" },          // infinite bounds
        { "1/(x^4+1)", "0", "inf", "1/4 * sqrt(2) * pi" },
        { "x/(x^3+1)", "0", "inf", "2/9 * sqrt(3) * pi" },               // ln terms: lc ratio 1
        { "1/(x^6+2)", "-inf", "inf", "1/3 * 2^(1/6) * pi" },
        { "1/(x^2+2x+5)", "-inf", "1", "3/8 * pi" },
        { "1/(x^2+1)", "inf", "0", "-1/2 * pi" },
        { "1/((x-1)(x-2))", "3", "inf", "ln(2)" },                       // ln|x-2| - ln|x-1| → 0
        { "1/(4x^2-1)", "1", "inf", "1/4 * ln(3)" },                     // lc ratio ln(2x-1)/(2x+1)
        { "1/(x^2+1)^2", "0", "inf", "1/4 * pi" },                       // rational term → 0
        { "x/(x^4+1)", "-inf", "0", "-1/4 * pi" },                       // atan(x^2) → +π/2 at -∞
        { "(x^2+1)/(x^4+3x^2+1)", "-inf", "inf", "2/5 * sqrt(5) * pi" }, // atan jump over R
        { "1/(x^2+1)", "inf", "inf", "0" },
        { "1/((x-1000)^2+1)", "-inf", "inf", "pi" },                    // peak far from 0: cut there
        { "1/((x-50)^2+1/10000)", "0", "100", "200 * atan(5000)" },     // narrow peak, no cancellation
        { "-8/(9x^2+12x+9)", "-2/3", "inf", "-4/15 * sqrt(5) * pi" },  // critical point at a bound
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        DefiniteResult r = definite_str(cases[i][0], cases[i][1], cases[i][2]);
        assert(!r.error && r.exact);
        expect_str(expr_to_string(r.exact), cases[i][3]);
        definite_free(&r);
    }
    DefiniteResult d = definite_str("1/x", "-1", "1");               // pole inside
    assert(d.error && !d.exact);
    d = definite_str("1/(x-1)^2", "1", "2");                         // pole at a bound
    assert(d.error && !d.exact);
    d = definite_str("1/(x^5+x+3)", "0", "1");                       // RootSum: number only
    assert(!d.error && !d.exact && d.note && fabs(d.approx - 0.2777409846) < 1e-9);
    d = definite_str("x/(x^2+1)", "0", "inf");                       // decays only like 1/x
    assert(d.error && strstr(d.error, "infinity"));
    d = definite_str("1/(x^2-1)", "2", "-inf");                      // pole on the way to -∞
    assert(d.error && strstr(d.error, "pole"));
    d = definite_str("x^2", "inf", "inf");                           // still diverges, not 0
    assert(d.error);
    d = definite_str("x^2/(x^4+x^2+3)", "-inf", "inf");              // RootSum: number only
    assert(!d.error && !d.exact && fabs(d.approx - 1.486903682) < 1e-8);
    printf("  test_reg_integrate_definite PASSED (%zu exact, 7 special)\n", sizeof cases / sizeof cases[0]);
}

// expr_lead: the leading coefficient and degree of the expression shapes that appear in
// antiderivatives (sums with the higher degree on either side, surd products, quotients,
// x-free powers), and NULL when the leading terms cancel
static Expr* sqrt_of(long m) { return create_sqrt(create_const(m)); }

static void check_lead(Expr* e, const char* coef, int deg) {
    int d;
    Expr* c = expr_lead(e, &d);
    if (!coef) assert(!c);
    else { assert(c && d == deg); expect_str(expr_to_string(c), coef); free_expr(c); }
    free_expr(e);
}

static void test_reg_expr_lead(void) {
    Expr* x = create_var("x");
    check_lead(create_sub(create_const(1), create_mul(sqrt_of(2), expr_copy(x))), "-sqrt(2)", 1);
    check_lead(create_add(create_pow(expr_copy(x), create_const(2)), create_mul(sqrt_of(3), expr_copy(x))), "1", 2);
    check_lead(create_div(create_sub(create_mul(create_const(2), expr_copy(x)), create_const(1)), sqrt_of(3)),
               "2/3 * sqrt(3)", 1);
    check_lead(create_div(create_const(5), create_add(expr_copy(x), create_const(1))), "5", -1);
    check_lead(create_mul(create_pow(create_const(2), create_const_str("0.5")), expr_copy(x)), "2^(1/2)", 1);
    check_lead(create_neg(create_pow(create_sub(create_const(3), expr_copy(x)), create_const(3))), "1", 3);
    check_lead(create_sub(create_mul(sqrt_of(2), expr_copy(x)), create_mul(sqrt_of(2), expr_copy(x))), NULL, 0);
    free_expr(x);
    printf("  test_reg_expr_lead PASSED (7 cases)\n");
}

// Q(sqrt(m)) arithmetic: parse, the four operations, exact signs at near-ties, printing
static void qs_of(QSurd* a, long p, long q, long m) {
    Expr* e = create_add(create_const(p), create_mul(create_const(q), sqrt_of(m)));
    assert(qs_parse(e, a));
    free_expr(e);
}

static char* expr_to_string_free(Expr* e) {
    char* s = expr_to_string(e);
    free_expr(e);
    return s;
}

static void expect_qs(const QSurd* a, const char* s) {
    Expr* e = qs_expr(a);
    expect_str(expr_to_string(e), s);
    free_expr(e);
}

// square roots in antiderivative constants merge (sqrt(2)·sqrt(3) = sqrt(6)); beside a radical
// that is not a square root they stay apart (2^(1/12)·sqrt(3)); trig values stay symbolic
static void test_reg_radical_products(void) {
    const char* cases[][2] = {
        { "1/(x^6+8)", "1/24 * sqrt(2) * atan(x / sqrt(2)) - 1/96 * sqrt(6) * ln(abs(x^2 - sqrt(6) * x + 2)) + 1/48 * sqrt(2) * atan(sqrt(2) * x - sqrt(3)) + 1/96 * sqrt(6) * ln(abs(x^2 + sqrt(6) * x + 2)) + 1/48 * sqrt(2) * atan(sqrt(2) * x + sqrt(3))" },
        { "1/(x^6+216)", "1/648 * sqrt(6) * atan(x / sqrt(6)) - 1/864 * sqrt(2) * ln(abs(x^2 - 3 * sqrt(2) * x + 6)) + 1/1296 * sqrt(6) * atan(1/3 * sqrt(6) * x - sqrt(3)) + 1/864 * sqrt(2) * ln(abs(x^2 + 3 * sqrt(2) * x + 6)) + 1/1296 * sqrt(6) * atan(1/3 * sqrt(6) * x + sqrt(3))" },                     // sqrt(6)·sqrt(3) = 3·sqrt(2)
        { "1/(x^12-2)", "1/24 * 2^(1/12) * ln(abs(x - 2^(1/12))) + 1/48 * 2^(1/12) * sqrt(3) * ln(abs(x^2 - 2^(1/12) * sqrt(3) * x + 2^(1/6))) - 1/24 * 2^(1/12) * atan(2^(11/12) * x - sqrt(3)) + 1/48 * 2^(1/12) * ln(abs(x^2 - 2^(1/12) * x + 2^(1/6))) - 1/24 * 2^(1/12) * sqrt(3) * atan(1/3 * 2^(11/12) * sqrt(3) * x - 1/3 * sqrt(3)) - 1/12 * 2^(1/12) * atan(1/2 * 2^(11/12) * x) - 1/48 * 2^(1/12) * ln(abs(x^2 + 2^(1/12) * x + 2^(1/6))) - 1/24 * 2^(1/12) * sqrt(3) * atan(1/3 * 2^(11/12) * sqrt(3) * x + 1/3 * sqrt(3)) - 1/48 * 2^(1/12) * sqrt(3) * ln(abs(x^2 + 2^(1/12) * sqrt(3) * x + 2^(1/6))) - 1/24 * 2^(1/12) * atan(2^(11/12) * x + sqrt(3)) - 1/24 * 2^(1/12) * ln(abs(x + 2^(1/12)))" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        Expr* f = parse_and_simplify(cases[i][0]);
        Expr* F = integrate(f, "x");
        expect_str(expr_to_string(F), cases[i][1]);
        free_expr(f); free_expr(F);
    }
    Expr* f = parse_and_simplify("1/(x^5+1)");
    Expr* F = integrate(f, "x");
    char* str = expr_to_string(F);
    assert(strstr(str, "cos(pi/5)") && !strstr(str, "sqrt(5)"));    // kept as cos(pi/5)
    free(str); free_expr(f); free_expr(F);
    printf("  test_reg_radical_products PASSED\n");
}

static void test_reg_qsurd(void) {
    QSurd a, b, r;
    qs_init(&a); qs_init(&b); qs_init(&r);
    qs_of(&a, 2, 1, 2); qs_of(&b, 2, -1, 2);
    assert(qs_div(&r, &a, &b)); expect_qs(&r, "3 + 2 * sqrt(2)");     // (2+√2)/(2-√2)
    assert(qs_mul(&r, &a, &b)); expect_qs(&r, "2");                   // norm
    assert(qs_sub(&r, &a, &b)); expect_qs(&r, "2 * sqrt(2)");
    qs_of(&a, 1, -1, 2); assert(qs_sign(&a) < 0);                     // 1 - √2
    qs_of(&a, 7, -4, 3); assert(qs_sign(&a) > 0);                     // 49 > 48
    qs_of(&a, -7, 4, 3); assert(qs_sign(&a) < 0);
    qs_of(&b, 1, 1, 3);
    qs_of(&a, 1, 1, 2); assert(!qs_mul(&r, &a, &b));                  // different fields
    Expr* e = create_div(create_mul(create_const(3), sqrt_of(5)), create_const(4));
    assert(qs_parse(e, &a)); expect_qs(&a, "3/4 * sqrt(5)");
    free_expr(e);
    e = create_mul(sqrt_of(2), sqrt_of(3));                           // two radicals: not one field
    assert(!qs_parse(e, &a));
    free_expr(e);
    qs_of(&a, 3, 2, 2); assert(qs_sqrt(&r, &a)); expect_qs(&r, "1 + sqrt(2)");   // (1+√2)²
    qs_of(&a, 7, -4, 3); assert(qs_sqrt(&r, &a)); expect_qs(&r, "2 - sqrt(3)");  // positive root
    qs_of(&a, 2, 1, 2); assert(!qs_sqrt(&r, &a));
    const long atans[][3] = { { -1, 1, 2 }, { 1, 1, 2 }, { 2, -1, 3 }, { 2, 1, 3 }, { -2, -1, 3 } };
    const char* turns[] = { "1/8 * pi", "3/8 * pi", "1/12 * pi", "5/12 * pi", "-5/12 * pi" };
    Rational zero = rat_zero_val();
    for (int i = 0; i < 5; i++) {                                     // tan of π/8, 3π/8, π/12, 5π/12
        Expr* t = create_atan(create_add(create_const(atans[i][0]), create_mul(create_const(atans[i][1]), sqrt_of(atans[i][2]))));
        Expr* v = expr_at(t, &zero);
        expect_str(expr_to_string(v), turns[i]);
        free_expr(t); free_expr(v);
    }
    rat_clear(&zero);
    qs_of(&a, 9, -4, 2); assert(qs_sqrt(&r, &a)); expect_qs(&r, "-1 + 2 * sqrt(2)"); // not 1 - 2√2
    qs_clear(&a); qs_clear(&b); qs_clear(&r);
    Expr* sums[][2] = {
        { create_atan(create_const(-2)), create_atan(create_div(create_const(-1), create_const(2))) }, // 1 - uv = 0, u < 0
        { create_atan(create_const(2)), create_atan(create_const(3)) },            // 1 - uv < 0: + π
        { create_mul(sqrt_of(2), create_pi()), create_mul(sqrt_of(2), create_ln(create_const(3))) },
    };
    const char* folded[] = { "-1/2 * pi", "3/4 * pi", "sqrt(2) * ln(3) + sqrt(2) * pi" };
    for (int i = 0; i < 3; i++) expect_str(expr_to_string_free(const_fold_sum(sums[i], 2)), folded[i]);
    Expr* lone[][1] = {
        { create_ln(create_add(create_const(3), create_mul(create_const(2), sqrt_of(2)))) },   // (1+√2)²
        { create_ln(create_add(create_const(-1), create_mul(create_const(3), sqrt_of(2)))) },  // conjugate
    };
    const char* canon[] = { "2 * ln(1 + sqrt(2))", "-ln(1 + 3 * sqrt(2)) + ln(17)" };
    for (int i = 0; i < 2; i++) expect_str(expr_to_string_free(const_fold_sum(lone[i], 1)), canon[i]);
    Expr* ratio[] = {                                                 // c_j/c_i = -1/2: merged as j, i
        create_mul(create_div(sqrt_of(2), create_const(2)), create_ln(create_add(create_const(1), sqrt_of(2)))),
        create_mul(create_div(sqrt_of(2), create_const(-4)), create_ln(create_add(create_const(1), create_mul(create_const(3), sqrt_of(2))))),
    };
    expect_str(expr_to_string_free(const_fold_sum(ratio, 2)), "-1/4 * sqrt(2) * ln(17) + 1/4 * sqrt(2) * ln(9 + 7 * sqrt(2))");
    Expr* roots[][3] = {                                              // square roots multiply into one
        { create_mul(sqrt_of(2), sqrt_of(6)), create_const(0), create_const(0) },
        { create_mul(sqrt_of(3), sqrt_of(3)), create_const(0), create_const(0) },
        { create_mul(create_mul(sqrt_of(2), sqrt_of(3)), sqrt_of(5)), create_const(0), create_const(0) },
    };
    const char* merged[] = { "2 * sqrt(3)", "3", "sqrt(30)" };
    for (int i = 0; i < 3; i++) expect_str(expr_to_string_free(const_fold_sum(roots[i], 3)), merged[i]);
    Expr* rest[] = { create_mul(sqrt_of(2), create_pi()), sqrt_of(3) };   // π multiples last
    expect_str(expr_to_string_free(const_fold_sum(rest, 2)), "sqrt(3) + sqrt(2) * pi");
    e = create_ln(create_add(create_const(1), create_mul(create_const(-1), sqrt_of(2))));
    expect_str(expr_to_string(e), "ln(abs(1 - sqrt(2)))");            // negative: keeps abs
    free_expr(e);
    printf("  test_reg_qsurd PASSED\n");
}

// lrt_log_part on its own, including reducible square-free d where R(t) has repeated
// roots (x^3 - x with a = 1: coefficients -1, 1/2, 1/2 force S_2 = a - t·d')
static void test_lrt_direct(void) {
    const char* cases[][2] = {
        { "1", "x^3 - x" }, { "x + 2", "x^3 - x" }, { "1", "(x-1)(x-2)(x-3)" },
        { "5x^2 - 1", "x^3 - x" }, { "1", "x^4 - 5x^2 + 4" }, { "x^2", "(x^2+1)(x^2+4)" },
        { "x^3 + 1", "(x^3-2)(x^2+1)" }, { "1", "x^5 + x + 3" }, { "x^2 + x", "x^4 + 3" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        Expr *a = parse_and_simplify(cases[i][0]), *d = parse_and_simplify(cases[i][1]);
        assert(a && d && a->type == EXPR_POLY && d->type == EXPR_POLY);
        Expr* F = lrt_log_part(a->poly, d->poly);
        if (!F || !derivative_is(F, a->poly, d->poly)) {
            char* Fs = F ? expr_to_string(F) : NULL;
            fprintf(stderr, "    LRT failed for %s / (%s): %s\n", cases[i][0], cases[i][1], Fs ? Fs : "NULL");
            free(Fs);
            assert(0);
        }
        free_expr(F); free_expr(a); free_expr(d);
    }
    printf("  test_lrt_direct PASSED (%zu cases)\n", sizeof cases / sizeof cases[0]);
}

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

// isolation may hand back an interval whose endpoint is a neighbouring exact root, e.g.
// (0, 4) beside [0, 0] for x(x-1)(x+1); refinement must still close in on the inner root
static void test_reg_refine_shared_endpoint(void) {
    const char* polys[] = { "x(x-1)(x+1)", "(x-1)(x-2)(x^2-3)", "x(x^2-2)", "(x-1)(x^3-2)" };
    Rational w = rat_from_si(1, 1000000), vlo, vhi;
    rat_init(&vlo); rat_init(&vhi);
    for (size_t i = 0; i < sizeof polys / sizeof polys[0]; i++) {
        Expr* e = parse_and_simplify(polys[i]);
        RootIntervalList L = isolate_real_roots(e->poly);
        for (int k = 0; k < L.count; k++) {
            refine_root(&L.items[k], e->poly, &w);
            poly_eval_rat(&vlo, e->poly, &L.items[k].lo);
            poly_eval_rat(&vhi, e->poly, &L.items[k].hi);
            bool exact = rat_eq(&L.items[k].lo, &L.items[k].hi) && rat_is_zero(&vlo);
            bool narrow = rat_to_double(&L.items[k].hi) - rat_to_double(&L.items[k].lo) <= 1e-6;
            assert(exact || (narrow && mpq_sgn(vlo.val) * mpq_sgn(vhi.val) < 0));
        }
        root_interval_list_free(&L);
        free_expr(e);
    }
    rat_clear(&w); rat_clear(&vlo); rat_clear(&vhi);
    printf("  test_reg_refine_shared_endpoint PASSED\n");
}

static void test_reg_solve_irrational_exact(void) {
    Expr* e = parse_and_simplify("x^3 - 2");
    Solutions* s = solve_polynomial(e->poly);
    assert(s->num_rational == 0 && s->num_irr == 1 && s->irr[0].num_real == 1);
    double lo = rat_to_double(&s->irr[0].real[0].lo), hi = rat_to_double(&s->irr[0].real[0].hi);
    assert(lo <= cbrt(2.0) && cbrt(2.0) <= hi && hi - lo < 1e-11);
    char* str = solutions_to_string(s);
    assert(strstr(str, "x = 2^(1/3) ~ 1.2599210498"));                    // exact, plus the decimal
    assert(strstr(str, "x = -1/2 * 2^(1/3) +/- 1/2 * 2^(1/3) * sqrt(3) * i"));
    free(str);
    free_solutions(s);
    free_expr(e);
    Expr* q = parse_and_simplify("x^2 + x + 1");
    s = solve_polynomial(q->poly);
    expect_str(solutions_to_string(s), "x = -1/2 +/- 1/2*sqrt(3)*i");
    free_solutions(s);
    free_expr(q);
    const char* known[][2] = {                                              // binomial, cyclotomic
        { "x^6 + 2", "x = 1/2 * 2^(1/6) * sqrt(3) +/- 1/2 * 2^(1/6) * i ~ (0.9720806486 +/- 0.5612310242i), "
                     "x = +/- 2^(1/6) * i ~ (0.0000000000 +/- 1.1224620483i), "
                     "x = -1/2 * 2^(1/6) * sqrt(3) +/- 1/2 * 2^(1/6) * i ~ (-0.9720806486 +/- 0.5612310242i)" },
        { "x^5 - 1", "x = 1, x = cos(2*pi/5) +/- sin(2*pi/5) * i ~ (0.3090169944 +/- 0.9510565163i), "
                     "x = -cos(pi/5) +/- sin(pi/5) * i ~ (-0.8090169944 +/- 0.5877852523i)" },
    };
    for (size_t i = 0; i < sizeof known / sizeof known[0]; i++) {
        q = parse_and_simplify(known[i][0]);
        s = solve_polynomial(q->poly);
        expect_str(solutions_to_string(s), known[i][1]);
        free_solutions(s);
        free_expr(q);
    }
    printf("  test_reg_solve_irrational_exact PASSED\n");
}

static void test_reg_integrate_complete(void) {
    Expr* r = integrate_str("1/((x-1)(x^2+1))"); // used to silently drop two terms
    expect_str(expr_to_string(r), "1/2 * ln(abs(x - 1)) - 1/4 * ln(abs(x^2 + 1)) - 1/2 * atan(x)");
    free_expr(r);
    r = integrate_str("1/(2x^2 - 3x + 1)");      // used to come out doubled
    expect_str(expr_to_string(r), "-ln(abs(2x - 1)) + ln(abs(x - 1))");
    free_expr(r);
    r = integrate_str("1/(x^2+1)^2");
    expect_str(expr_to_string(r), "1/2x / (x^2 + 1) + 1/2 * atan(x)");
    free_expr(r);
    printf("  test_reg_integrate_complete PASSED\n");
}

/***********************************************************************
*  property tests on generated inputs
*************************************************************************/

static long gcd_l(long a, long b) { return b == 0 ? a : gcd_l(b, a % b); }

// a random factor whose irreducibility over Q is known by construction: linear a·x + b,
// x^2 + k (k > 0), x^2 - d (d not a square), a rootless cubic, or x^4 + k; returns real-root count
// x^3 + b·x + c with no integer root dividing c is irreducible over Q
static Polynomial* irreducible_cubic(int* real_roots) {
    long c[4] = { 0, 0, 0, 1 };
    bool has_root;
    do {
        c[1] = rnd(-6, 6);
        c[0] = rnd(-9, 9);
        has_root = c[0] == 0;
        for (long r = 1; r <= labs(c[0]) && !has_root; r++)
            if (c[0] % r == 0)
                for (long sgn = -1; sgn <= 1; sgn += 2) {
                    long x = sgn * r;
                    if (x * x * x + c[1] * x + c[0] == 0) has_root = true;
                }
    } while (has_root);
    long disc = -4 * c[1] * c[1] * c[1] - 27 * c[0] * c[0]; // > 0: three real roots
    *real_roots = disc > 0 ? 3 : 1;
    return poly_si(c, 3);
}

static Polynomial* known_irreducible(int* real_roots) {
    long c[3] = { 0, 0, 0 };
    int kind = (int)rnd(0, 4);
    if (kind == 3) return irreducible_cubic(real_roots);
    if (kind == 4) { // x^4 + k is irreducible for k > 0 unless k = 4c^4
        static const long ks[] = { 1, 2, 3, 5, 6, 7 };
        long q[5] = { ks[rnd(0, 5)], 0, 0, 0, 1 };
        *real_roots = 0;
        return poly_si(q, 4);
    }
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

// k·(U·W' - U'·W) / (U^2 + shift·W^2): the derivative of every atan / log-ratio shape
typedef struct { const Rational* k; const Rational* shift; } CrossSpec;

static RF rf_cross(const Polynomial* U, const Polynomial* W, const CrossSpec* sp) {
    Polynomial *up = poly_derivative(U), *wp = poly_derivative(W);
    Polynomial *a = poly_multiply(U, wp), *b = poly_multiply(up, W), *n = poly_subtract(a, b);
    Polynomial *u2 = poly_multiply(U, U), *w2 = poly_multiply(W, W), *sw2 = poly_scale(w2, sp->shift);
    RF r = rf_make(poly_scale(n, sp->k), poly_add(u2, sw2));
    Polynomial* tmp[] = { up, wp, a, b, n, u2, w2, sw2 };
    for (size_t i = 0; i < sizeof tmp / sizeof tmp[0]; i++) free_polynomial(tmp[i]);
    return r;
}

// "sqrt(m)" (W = 1) or "sqrt(m) * W": returns W (caller frees) or NULL for another shape
static Polynomial* sqrt_w(const Expr* e) {
    if (sqrt_arg(e)) return poly_one_t();
    if (e->type == EXPR_BINOP && e->binop.op == OP_MUL && sqrt_arg(e->binop.left) && e->binop.right->type == EXPR_POLY)
        return poly_copy(e->binop.right->poly);
    return NULL;
}

// sqrt(m)·h for h = atan(U/(√m·W)), atan(√m·W/U) or ln|(U ± √m·W)/(U ∓ √m·W)|
static bool rf_sqrt_term(const Expr* rest, const Rational* m, RF* out) {
    Rational k, sh;
    rat_init(&k); rat_init(&sh);
    bool ok = false;
    const Expr* arg = rest->type == EXPR_ATAN ? rest->fn_arg : rest->type == EXPR_LN ? rest->ln_arg : NULL;
    if (arg && arg->type == EXPR_BINOP && arg->binop.op == OP_DIV) {
        const Expr *l = arg->binop.left, *r = arg->binop.right;
        Polynomial* Wr = (rest->type == EXPR_ATAN && l->type == EXPR_POLY) ? sqrt_w(r) : NULL;
        if (Wr) {                                                            // atan(U/(√m·W)):
            rat_neg(&k, m); rat_set(&sh, m);                                 // -m·(UW' - U'W)/(U^2 + mW^2)
            CrossSpec sp = { &k, &sh };
            *out = rf_cross(l->poly, Wr, &sp);
            free_polynomial(Wr);
            ok = true;
        } else if (rest->type == EXPR_ATAN && r->type == EXPR_POLY) {        // m·(UW' - U'W)/(U^2 + mW^2)
            Polynomial* W = sqrt_w(l);
            if (W) {
                rat_set(&k, m); rat_set(&sh, m);
                CrossSpec sp = { &k, &sh };
                *out = rf_cross(r->poly, W, &sp);
                free_polynomial(W);
                ok = true;
            }
        } else if (rest->type == EXPR_LN && l->type == EXPR_BINOP && l->binop.left->type == EXPR_POLY) {
            Polynomial* W = sqrt_w(l->binop.right);                          // ±2m·(UW' - U'W)/(U^2 - mW^2)
            if (W) {
                rat_add(&k, m, m);
                if (l->binop.op == OP_SUB) rat_neg(&k, &k);
                rat_neg(&sh, m);
                CrossSpec sp = { &k, &sh };
                *out = rf_cross(l->binop.left->poly, W, &sp);
                free_polynomial(W);
                ok = true;
            }
        }
    }
    rat_clear(&k); rat_clear(&sh);
    return ok;
}

// derivative of one emitted term c·(...)
static bool rf_scaled_term(const Expr* t, RF* out) {
    const Rational* c = &t->binop.left->constant;
    const Expr* g = t->binop.right;
    bool ok = true;
    Rational one = rat_one_val();
    if (g->type == EXPR_LN && g->ln_arg->type == EXPR_POLY) {
        *out = rf_make(poly_derivative(g->ln_arg->poly), poly_copy(g->ln_arg->poly));
    } else if (g->type == EXPR_ATAN && g->fn_arg->type == EXPR_POLY) {
        *out = rf_arc(g->fn_arg->poly, &one, &one);
    } else if (g->type == EXPR_ATAN && g->fn_arg->type == EXPR_BINOP && g->fn_arg->binop.op == OP_DIV &&
               g->fn_arg->binop.left->type == EXPR_POLY && g->fn_arg->binop.right->type == EXPR_POLY) {
        CrossSpec sp = { &one, &one };                                       // atan(N/G)
        *out = rf_cross(g->fn_arg->binop.right->poly, g->fn_arg->binop.left->poly, &sp);
    } else if (g->type == EXPR_BINOP && g->binop.op == OP_DIV) {
        *out = rf_quotient_rule(g->binop.left->poly, g->binop.right->poly);
    } else if (g->type == EXPR_BINOP && g->binop.op == OP_MUL && sqrt_arg(g->binop.left)) {
        ok = rf_sqrt_term(g->binop.right, sqrt_arg(g->binop.left), out);
    } else {
        ok = false;
    }
    rat_clear(&one);
    if (ok) rf_scale(out, c);
    return ok;
}

typedef struct {
    const Expr* items[16];
    int count;
} RootSums;

// exact derivative of F into acc; RootSum terms are set aside for a numeric check
static bool rf_derivative(const Expr* F, RF* acc, RootSums* rs) {
    RF t;
    if (F->type == EXPR_BINOP && F->binop.op == OP_ADD)
        return rf_derivative(F->binop.left, acc, rs) && rf_derivative(F->binop.right, acc, rs);
    if (F->type == EXPR_ROOTSUM && rs->count < 16) { rs->items[rs->count++] = F; return true; }
    if (F->type == EXPR_POLY) t = rf_make(poly_derivative(F->poly), poly_one_t());
    else if (F->type == EXPR_BINOP && F->binop.op == OP_DIV) t = rf_quotient_rule(F->binop.left->poly, F->binop.right->poly);
    else if (F->type == EXPR_BINOP && F->binop.op == OP_MUL && F->binop.left->type == EXPR_CONST) {
        if (!rf_scaled_term(F, &t)) return false;
    } else return false;
    rf_add_into(acc, t);
    return true;
}

// --- exact RootSum derivative at a rational point, via traces in K = Q[t]/(q) ---

// power sums p_i = Σ alpha^i over the roots of q, i < deg q (Newton's identities)
static Rational* power_sums(const Polynomial* q) {
    int n = polynomial_degree(q);
    Rational *a = malloc((size_t)(n + 1) * sizeof(Rational)), *ps = malloc((size_t)n * sizeof(Rational));
    Rational lc, t;
    rat_init(&lc); rat_init(&t);
    poly_coeff(&lc, q, n);
    for (int k = 0; k <= n; k++) { rat_init(&a[k]); poly_coeff(&a[k], q, k); rat_div(&a[k], &a[k], &lc); }
    for (int k = 0; k < n; k++) {
        rat_init(&ps[k]);
        if (k == 0) { mpq_set_si(ps[0].val, n, 1); continue; }
        Rational kk = rat_from_int(k);                // p_k = -(k·a_{n-k} + Σ a_{n-j}·p_{k-j})
        rat_mul(&ps[k], &kk, &a[n - k]);
        for (int j = 1; j < k; j++) { rat_mul(&t, &a[n - j], &ps[k - j]); rat_add(&ps[k], &ps[k], &t); }
        rat_neg(&ps[k], &ps[k]);
        rat_clear(&kk);
    }
    for (int k = 0; k <= n; k++) rat_clear(&a[k]);
    free(a);
    rat_clear(&lc); rat_clear(&t);
    return ps;
}

// S(t, x0) (deriv = 0) or dS/dx(t, x0) (deriv = 1) as a polynomial in t
static Polynomial* s_at(const Expr* e, const Rational* x0, int deriv) {
    Polynomial* r = create_polynomial();
    Rational xp, w;
    rat_init(&xp); rat_init(&w);
    for (int k = deriv; k <= e->rootsum.s_deg; k++) {
        rat_pow_int(&xp, x0, k - deriv);
        if (deriv) { Rational kk = rat_from_int(k); rat_mul(&xp, &xp, &kk); rat_clear(&kk); }
        Polynomial* term = poly_scale(e->rootsum.s[k], &xp);
        Polynomial* sum2 = poly_add(r, term);
        free_polynomial(r); free_polynomial(term);
        r = sum2;
    }
    rat_clear(&xp); rat_clear(&w);
    return r;
}

// Σ alpha·S_x(alpha, x0)/S(alpha, x0) = Tr(t·S_x·S^-1 mod q), exactly; false at a pole
static bool rootsum_derivative_at(const Expr* e, const Rational* x0, Rational* out) {
    const Polynomial* q = e->rootsum.q;
    Polynomial *sv = s_at(e, x0, 0), *sx = s_at(e, x0, 1);
    PolyBezout bz = poly_xgcd(sv, q);
    bool ok = polynomial_degree(bz.g) == 0;           // S(alpha, x0) != 0 for every root
    if (ok) {
        Rational one = rat_one_val();
        Polynomial* tpoly = create_polynomial();
        add_term(tpoly, &one, 1);
        Polynomial *u = poly_multiply(tpoly, sx), *v = poly_multiply(u, bz.s), *el = poly_mod(v, q);
        Rational* ps = power_sums(q);
        Rational c;
        rat_init(&c);
        mpq_set_si(out->val, 0, 1);
        for (int i = 0; i < polynomial_degree(q); i++) {
            poly_coeff(&c, el, i);
            rat_mul(&c, &c, &ps[i]);
            rat_add(out, out, &c);
            rat_clear(&ps[i]);
        }
        free(ps);
        rat_clear(&c); rat_clear(&one);
        free_polynomial(tpoly); free_polynomial(u); free_polynomial(v); free_polynomial(el);
    }
    free_bezout(&bz);
    free_polynomial(sv); free_polynomial(sx);
    return ok;
}

// F' == num/den: as rational functions when F has no RootSum, else exactly at rational points
static bool derivative_is(const Expr* F, const Polynomial* num, const Polynomial* den) {
    RF acc = rf_make(create_polynomial(), poly_one_t());
    RootSums rs = { { NULL }, 0 };
    bool ok = rf_derivative(F, &acc, &rs);
    if (ok && rs.count == 0) {
        Polynomial *x = poly_multiply(acc.num, den), *y = poly_multiply(num, acc.den);
        ok = poly_equal(x, y);
        free_polynomial(x); free_polynomial(y);
    }
    static const long pts[][2] = { { -10, 3 }, { -17, 10 }, { -3, 5 }, { 9, 20 }, { 13, 10 }, { 29, 10 } };
    Rational x0, f, d, got, part;
    rat_init(&x0); rat_init(&f); rat_init(&d); rat_init(&got); rat_init(&part);
    int checked = 0;
    for (size_t i = 0; ok && rs.count > 0 && i < sizeof pts / sizeof pts[0]; i++) {
        rat_set_si(&x0, pts[i][0], pts[i][1]);
        poly_eval_rat(&d, den, &x0);
        poly_eval_rat(&part, acc.den, &x0);
        if (rat_is_zero(&d) || rat_is_zero(&part)) continue; // a pole: skip this point
        poly_eval_rat(&f, num, &x0); rat_div(&f, &f, &d);
        poly_eval_rat(&got, acc.num, &x0); rat_div(&got, &got, &part);
        bool pole = false;
        for (int k = 0; k < rs.count && !pole; k++) {
            pole = !rootsum_derivative_at(rs.items[k], &x0, &part);
            rat_add(&got, &got, &part);
        }
        if (pole) continue;
        ok = rat_eq(&got, &f);
        checked++;
    }
    if (rs.count > 0 && checked < 3) ok = false;        // demand real evidence
    rat_clear(&x0); rat_clear(&f); rat_clear(&d); rat_clear(&got); rat_clear(&part);
    rf_free(&acc);
    return ok;
}

// integrate()'s user-facing answer: numerically checked when RootSum was expanded; true if so
static bool check_expanded(Expr* integrand, const Polynomial* num, const Polynomial* den) {
    Expr* G = integrate(integrand, "x");
    char* Gs = expr_to_string(G);
    bool radical = !strstr(Gs, "RootSum") && (strstr(Gs, "^(") || strstr(Gs, "pi"));
    if (radical && !numeric_derivative_is(G, num, den)) {
        fprintf(stderr, "    expanded antiderivative wrong: %s\n", Gs);
        assert(0);
    }
    free(Gs);
    free_expr(G);
    return radical;
}

// backward generation: integrands built from known irreducible factors (incl. cubics -> LRT)
static void test_prop_integrate_derivative(void) {
    int with_rootsum = 0, expanded = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        Generated g = generate_product();
        long nc[6];
        int nd = (int)rnd(0, polynomial_degree(g.product) + 1);
        if (nd > 5) nd = 5;
        for (int i = 0; i <= nd; i++) nc[i] = rnd(-9, 9);
        nc[nd] = nc[nd] ? nc[nd] : 1;
        Polynomial* num = poly_si(nc, nd);
        Expr* integrand = create_div(create_poly_expr(poly_copy(num)), create_poly_expr(poly_copy(g.product)));
        Expr* F = integrate_rootsum(integrand, "x");      // exact form, checked exactly below
        if (!F) fprintf(stderr, "    integrate failed: %s\n", integrate_error());
        assert(F != NULL);
        char* Fs0 = expr_to_string(F);
        with_rootsum += strstr(Fs0, "RootSum") != NULL;
        free(Fs0);
        if (!derivative_is(F, num, g.product)) {
            char *fs = expr_to_string(integrand), *Fs = expr_to_string(F);
            fprintf(stderr, "    d/dx F != f\n    f = %s\n    F = %s\n", fs, Fs);
            free(fs); free(Fs);
            assert(0);
        }
        expanded += check_expanded(integrand, num, g.product);
        free_expr(F);
        free_expr(integrand);
        free_polynomial(num);
        free_generated(&g);
    }
    printf("  test_prop_integrate_derivative PASSED (%d trials, %d with RootSum terms, %d expanded)\n",
           TRIALS, with_rootsum, expanded);
}

// random (a, d) with deg d in 2..7; when step > 1 both use only powers r + step·k, a
// symmetry that makes the chain over Q[t] defective (gaps between nonzero subresultants)
static RTPair random_rt_pair(void) {
    long dc[8] = { 0 }, ac[8] = { 0 };
    int step = (int)rnd(1, 3), m = step * (int)rnd(step == 1 ? 2 : 1, 7 / step), r = (int)rnd(0, step - 1);
    for (int i = 0; i < m; i += step) dc[i] = rnd(-3, 3);
    dc[m] = rnd(1, 3);
    for (int i = r; i < m; i += step) ac[i] = rnd(-3, 3);
    return (RTPair){ poly_si(dc, m), poly_si(ac, m - 1) };
}

// a == b as bivariate polynomials (missing coefficients are zero)
static bool ptx_equal(const PolyTX* a, const PolyTX* b) {
    int n = a->deg > b->deg ? a->deg : b->deg;
    for (int k = 0; k <= n; k++) {
        Polynomial* z = create_polynomial();
        bool eq = poly_equal(k <= a->deg ? a->c[k] : z, k <= b->deg ? b->c[k] : z);
        free_polynomial(z);
        if (!eq) return false;
    }
    return true;
}

// rt_chain(p) == the Sylvester-determinant definition, index by index; true if defective
static bool check_chain(const RTPair* p) {
    SubresChain ch = rt_chain(p);
    bool defective = false;
    for (int j = 0; j <= ch.m - 2; j++) {
        PolyTX ref = ref_subresultant(p, j);
        if (!ptx_equal(&ch.S[j], &ref)) {
            char *d = poly_to_string(p->d), *a = poly_to_string(p->a);
            fprintf(stderr, "    S_%d differs for d = %s, a = %s\n", j, d, a);
            assert(0);
        }
        defective |= !ptx_is_zero(&ch.S[j]) && ch.S[j].deg < j; // S_j similar to a lower one
        ptx_free(&ref);
    }
    chain_free(&ch);
    return defective;
}

// fixed chains with wide gaps (e.g. 1/(x^6+2): S_3 = S_2 = 0), then random pairs
static void test_prop_subres_chain(void) {
    const char* fixed[][2] = {
        { "1", "x^6 + 2" }, { "1", "x^4 + 1" }, { "x", "x^4 + 1" }, { "x^2", "x^6 + 1" },
        { "x^3", "x^8 + 5" }, { "1", "x^7 - 3" }, { "x + 1", "x^6 - x^3 + 2" }, { "0", "x^5 + x + 3" },
    };
    for (size_t i = 0; i < sizeof fixed / sizeof fixed[0]; i++) {
        Expr *a = parse_and_simplify(fixed[i][0]), *d = parse_and_simplify(fixed[i][1]);
        RTPair p = { d->poly, a->poly };
        check_chain(&p);
        free_expr(a); free_expr(d);
    }
    int defective = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        RTPair p = random_rt_pair();
        defective += check_chain(&p);
        free_polynomial((Polynomial*)p.d); free_polynomial((Polynomial*)p.a);
    }
    printf("  test_prop_subres_chain PASSED (%zu fixed, %d trials, %d defective)\n",
           sizeof fixed / sizeof fixed[0], TRIALS, defective);
}

// c·x^n ± e or x^(sN) ± b^N, sometimes times (x - k) or (x^2 + k)
static Polynomial* radical_denominator(void) {
    long dc[13] = { 0 }, sign = rnd(0, 1) ? 1 : -1;
    int n;
    if (rnd(0, 1)) {
        n = (int)rnd(3, 9);
        dc[n] = rnd(1, 4);
        dc[0] = sign * rnd(1, 7);
    } else {                                              // splits into rescaled cyclotomics
        int s = (int)rnd(1, 2), N = (int)rnd(3, 6);
        long b = rnd(1, 3), bN = 1;
        for (int i = 0; i < N; i++) bN *= b;
        n = s * N;
        dc[n] = 1;
        dc[0] = sign * bN;
    }
    Polynomial* d = poly_si(dc, n);
    long extra[3] = { rnd(-3, 3), rnd(0, 1), 0 };
    int shape = (int)rnd(0, 2);                           // 0: none, 1: x - k, 2: x^2 + k
    if (shape == 0) return d;
    if (shape == 1) extra[1] = 1;
    else { extra[0] = rnd(1, 5); extra[1] = 0; extra[2] = 1; }
    Polynomial* e = poly_si(extra, shape);
    Polynomial* p = poly_multiply(d, e);
    free_polynomial(d); free_polynomial(e);
    return p;
}

// denominators whose roots are known exactly: integrate() must expand RootSum, and the
// answer is checked numerically (exactly, when a RootSum is left)
static void test_prop_integrate_radicals(void) {
    int expanded = 0, kept = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        Polynomial* den = radical_denominator();
        long nc[8];
        int nd = (int)rnd(0, 6);
        for (int i = 0; i <= nd; i++) nc[i] = rnd(-5, 5);
        nc[nd] = nc[nd] ? nc[nd] : 1;
        Polynomial* num = poly_si(nc, nd);
        Expr* integrand = create_div(create_poly_expr(poly_copy(num)), create_poly_expr(poly_copy(den)));
        Expr* G = integrate(integrand, "x");
        assert(G != NULL);
        char* Gs = expr_to_string(G);
        bool has_rs = strstr(Gs, "RootSum") != NULL;
        bool ok = has_rs ? derivative_is(G, num, den) : numeric_derivative_is(G, num, den);
        if (!ok) {
            char* fs = expr_to_string(integrand);
            fprintf(stderr, "    d/dx F != f\n    f = %s\n    F = %s\n", fs, Gs);
            assert(0);
        }
        kept += has_rs;
        expanded += !has_rs && (strstr(Gs, "^(") || strstr(Gs, "pi"));
        free(Gs); free_expr(G); free_expr(integrand);
        free_polynomial(num); free_polynomial(den);
    }
    printf("  test_prop_integrate_radicals PASSED (%d trials, %d expanded, %d RootSum kept)\n",
           TRIALS, expanded, kept);
}

// P(z) in complex double precision
static double _Complex poly_at(const Polynomial* P, double _Complex z) {
    double _Complex v = 0;
    for (int k = polynomial_degree(P); k >= 0; k--) {
        Rational c;
        rat_init(&c);
        poly_coeff(&c, P, k);
        v = v * z + rat_to_double(&c);
        rat_clear(&c);
    }
    return v;
}

// every exact root that solve prints for a binomial/cyclotomic factor is a root of it,
// and there are exactly deg P of them (conjugate pairs counted twice)
static void test_prop_exact_roots(void) {
    int factors = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        Polynomial* d = radical_denominator();
        Solutions* s = solve_polynomial(d);
        for (int i = 0; i < s->num_irr; i++) {
            const Polynomial* P = s->irr[i].factor;
            RootSet rs;
            if (!known_roots(P, &rs)) continue;
            double scale = 0;
            for (const Term* t = P->head; t; t = t->next) scale += fabs(rat_to_double(&t->coeff));
            int roots = 0;
            for (int j = 0; j < rs.count; j++) {
                Rational theta = rat_from_si(rs.js[j], rs.L);
                Expr *re = root_part(&rs, &theta, TRIG_COS), *im = root_part(&rs, &theta, TRIG_SIN);
                bool ok = true;
                double a = re ? expr_eval(re, 0, &ok) : 0, b = im ? expr_eval(im, 0, &ok) : 0;
                assert(ok && cabs(poly_at(P, a + b * I)) < 1e-9 * scale);
                roots++;
                if (re) free_expr(re);
                if (im) free_expr(im);
                rat_clear(&theta);
            }
            assert(roots == polynomial_degree(P));
            root_set_free(&rs);
            factors++;
        }
        char* str = solutions_to_string(s);                   // and the printer copes with all of them
        assert(!strstr(str, "non-real roots of (") || strstr(str, "x^4"));  // only Capelli quartics left
        free(str);
        free_solutions(s);
        free_polynomial(d);
    }
    printf("  test_prop_exact_roots PASSED (%d trials, %d factors with exact roots)\n", TRIALS, factors);
}

// random rational bound p/q with |p| <= 12, q <= 4
static void random_bound(Rational* r) { rat_set_si(r, rnd(-12, 12), rnd(1, 4)); }

// definite integrals of generated integrands over random bounds: the exact value must always
// pass the built-in quadrature check (a failure leaves only the number and a note)
static void test_prop_integrate_definite(void) {
    int exact = 0, poles = 0, numeric = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        Generated g = generate_product();
        long nc[4];
        int nd = (int)rnd(0, 3);
        for (int i = 0; i <= nd; i++) nc[i] = rnd(-9, 9);
        nc[nd] = nc[nd] ? nc[nd] : 1;
        Polynomial* num = poly_si(nc, nd);
        Expr* f = create_div(create_poly_expr(poly_copy(num)), create_poly_expr(poly_copy(g.product)));
        Rational a, b;
        rat_init(&a); rat_init(&b);
        random_bound(&a);
        random_bound(&b);
        DefiniteResult r = integrate_definite(f, &(Bounds){ { &a, 0 }, { &b, 0 } });
        if (r.note && strstr(r.note, "internal")) {
            char* fs = expr_to_string(f);
            fprintf(stderr, "    definite integral check failed: %s on [%g, %g]\n", fs, rat_to_double(&a), rat_to_double(&b));
            assert(0);
        }
        exact += r.exact != NULL;
        poles += r.error != NULL;
        numeric += !r.error && !r.exact;
        definite_free(&r);
        rat_clear(&a); rat_clear(&b);
        free_expr(f); free_polynomial(num); free_generated(&g);
    }
    printf("  test_prop_integrate_definite PASSED (%d trials: %d exact, %d RootSum numeric, %d divergent)\n",
           TRIALS, exact, numeric, poles);
}

// p² + c with c > 0: a denominator with no real roots, so an infinite range has no pole
static Polynomial* positive_denominator(void) {
    long pc[4];
    int dp = (int)rnd(1, 3);
    for (int i = 0; i <= dp; i++) pc[i] = rnd(-4, 4);
    pc[dp] = pc[dp] ? pc[dp] : 1;
    Polynomial *p = poly_si(pc, dp), *sq = poly_multiply(p, p), *c = poly_si((long[]){ rnd(1, 6) }, 0);
    Polynomial* d = poly_add(sq, c);
    free_polynomial(p); free_polynomial(sq); free_polynomial(c);
    return d;
}

// infinite bounds on convergent integrands (deg num <= deg den - 2, no real pole) over
// (-∞, ∞), [a, ∞) and (-∞, b]: the exact limit must always pass the quadrature check
static void test_prop_integrate_infinite(void) {
    int exact = 0, numeric = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        Polynomial* den = positive_denominator();
        long nc[5];
        int nd = (int)rnd(0, polynomial_degree(den) - 2);
        for (int i = 0; i <= nd; i++) nc[i] = rnd(-9, 9);
        Polynomial* num = poly_si(nc, nd);
        Expr* f = create_div(create_poly_expr(num), create_poly_expr(den));
        Rational a, b;
        rat_init(&a); rat_init(&b);
        random_bound(&a);
        random_bound(&b);
        int kind = trial % 3;                                 // both infinite, then each side alone
        Bounds ab = { { &a, kind == 2 ? 0 : -1 }, { &b, kind == 1 ? 0 : 1 } };
        DefiniteResult r = integrate_definite(f, &ab);
        if (r.error || (r.note && strstr(r.note, "internal"))) {
            char* fs = expr_to_string(f);
            fprintf(stderr, "    infinite integral failed: %s (kind %d): %s\n", fs, kind, r.error ? r.error : r.note);
            assert(0);
        }
        assert(r.exact || strstr(r.note, "RootSum"));      // the limits themselves are always found
        exact += r.exact != NULL;
        numeric += !r.exact;
        definite_free(&r);
        rat_clear(&a); rat_clear(&b);
        free_expr(f);
    }
    printf("  test_prop_integrate_infinite PASSED (%d trials: %d exact, %d RootSum numeric)\n", TRIALS, exact, numeric);
}

// canonical definite values: an even integrand p(x²)/q(x²) prints the same on [a, b] and
// [-b, -a] (the antiderivative's terms arrive in another order with other arguments)
static void test_prop_definite_symmetric(void) {
    int compared = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        int k = (int)rnd(1, 4);
        long dc[9] = { 0 }, nc[9] = { 0 };
        dc[2 * k] = 1;
        dc[0] = rnd(-8, 8);
        if (rnd(0, 1)) dc[k] = rnd(-3, 3);                    // x^2k + b·x^k + c with k even below
        if (k % 2) dc[k] = 0;
        nc[0] = rnd(1, 4);
        if (k > 1 && rnd(0, 1)) nc[2] = rnd(-3, 3);
        Polynomial *num = poly_si(nc, k > 1 ? 2 : 0), *den = poly_si(dc, 2 * k);
        Expr* f = create_div(create_poly_expr(num), create_poly_expr(den));
        Rational a, b, na, nb;
        rat_init(&a); rat_init(&b); rat_init(&na); rat_init(&nb);
        rat_set_si(&a, rnd(0, 6), rnd(1, 3));
        rat_set_si(&b, rnd(7, 12), rnd(1, 3));
        rat_neg(&na, &a); rat_neg(&nb, &b);
        DefiniteResult r1 = integrate_definite(f, &(Bounds){ { &a, 0 }, { &b, 0 } });
        DefiniteResult r2 = integrate_definite(f, &(Bounds){ { &nb, 0 }, { &na, 0 } });
        if (r1.exact && r2.exact) {
            char *s1 = expr_to_string(r1.exact), *s2 = expr_to_string(r2.exact);
            if (strcmp(s1, s2) != 0) { fprintf(stderr, "    asymmetric:\n    %s\n    %s\n", s1, s2); assert(0); }
            free(s1); free(s2);
            compared++;
        }
        definite_free(&r1); definite_free(&r2);
        rat_clear(&a); rat_clear(&b); rat_clear(&na); rat_clear(&nb);
        free_expr(f);
    }
    printf("  test_prop_definite_symmetric PASSED (%d trials, %d exact pairs equal as strings)\n", TRIALS, compared);
}

// random polynomial of degree d, coefficients up to 10^big (or small rationals)
static Polynomial* random_q_poly(int d, int big) {
    Polynomial* p = create_polynomial();
    for (int i = 0; i <= d; i++) {
        Rational c;
        rat_init(&c);
        if (big) {
            char lit[64];
            snprintf(lit, sizeof lit, "%ld%09ld", rnd(1, 999999), rnd(0, 999999999));
            rat_parse_literal(&c, lit);
            if (rnd(0, 1)) rat_neg(&c, &c);
        } else rat_set_si(&c, rnd(-9, 9), rnd(1, 4));
        if (i == d && rat_is_zero(&c)) rat_set_si(&c, 1, 1);
        add_term(p, &c, i);
        rat_clear(&c);
    }
    return p;
}

static bool poly_eq(const Polynomial* a, const Polynomial* b) {
    Polynomial* d = poly_subtract(a, b);
    bool e = poly_is_zero(d);
    free_polynomial(d);
    return e;
}

// modular gcd and inverse agree with Euclid over Q: random a = g·u, b = g·v (so the gcd is
// usually g, with large or fractional coefficients), and inverses modulo random M
static void test_prop_modular_gcd(void) {
    int nontrivial = 0, inverses = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        int big = trial % 3 == 0;
        Polynomial *g = random_q_poly((int)rnd(0, 4), big), *u = random_q_poly((int)rnd(1, 6), big),
                   *v = random_q_poly((int)rnd(1, 6), big);
        Polynomial *a = poly_multiply(g, u), *b = poly_multiply(g, v);
        Polynomial *fast = modular_gcd(a, b), *slow = poly_gcd(a, b);
        assert(fast && poly_eq(fast, slow));
        nontrivial += polynomial_degree(slow) > 0;
        Polynomial* inv = modular_inverse(u, v);              // u invertible modulo v iff coprime
        if (polynomial_degree(v) > 0) {
            Polynomial* uv = poly_gcd(u, v);
            if (polynomial_degree(uv) == 0) {
                PolyBezout bz = poly_xgcd(u, v);
                Polynomial* s = poly_mod(bz.s, v);
                assert(inv && poly_eq(inv, s));
                free_polynomial(s);
                free_bezout(&bz);
                inverses++;
            }
            free_polynomial(uv);
        }
        free_polynomial(inv);
        free_polynomial(fast); free_polynomial(slow);
        free_polynomial(a); free_polynomial(b); free_polynomial(g); free_polynomial(u); free_polynomial(v);
    }
    // the first prime 2^31 - 1 is unlucky for x and x - (2^31 - 1) (both become x), and divides
    // the leading coefficient of (2^31 - 1)·x^2 + 1: the image is discarded or skipped
    Polynomial *x = poly_si((long[]){ 0, 1 }, 1), *y = poly_si((long[]){ -2147483647L, 1 }, 1);
    Polynomial* h = modular_gcd(x, y);
    assert(h && polynomial_degree(h) == 0);
    Polynomial* z = poly_si((long[]){ 1, 0, 2147483647L }, 2);
    Polynomial* hz = modular_gcd(z, x);
    assert(hz && polynomial_degree(hz) == 0);
    // x - p1·p2 looks like x for both first primes: a stable but wrong reconstruction that only
    // the exact check rejects
    Polynomial* y2 = poly_si((long[]){ -4611685975477714963L, 1 }, 1);
    Polynomial* h2 = modular_gcd(x, y2);
    assert(h2 && polynomial_degree(h2) == 0);
    // with 40-digit coefficients the gcd needs many primes, and the fifth (2147483563) is
    // unlucky for x and x - 2147483563: its longer image must be skipped, not mixed in
    Polynomial *g = random_q_poly(3, 1), *g2 = poly_multiply(g, g);
    Polynomial *gx = poly_multiply(g2, x), *y5 = poly_si((long[]){ -2147483563L, 1 }, 1), *gy = poly_multiply(g2, y5);
    Polynomial *h5 = modular_gcd(gx, gy), *m5 = poly_gcd(gx, gy);
    assert(h5 && poly_eq(h5, m5));
    free_polynomial(y2); free_polynomial(h2); free_polynomial(g); free_polynomial(g2); free_polynomial(gx);
    free_polynomial(y5); free_polynomial(gy); free_polynomial(h5); free_polynomial(m5);
    free_polynomial(x); free_polynomial(y); free_polynomial(h); free_polynomial(z); free_polynomial(hz);
    printf("  test_prop_modular_gcd PASSED (%d trials, %d nontrivial gcds, %d inverses)\n", TRIALS, nontrivial, inverses);
}

// --- exp integration (Risch) ---

static RischResult risch_str(const char* f) {
    Expr* e = parse_input(f);
    assert(e);
    RischResult r = integrate_exp(e);
    free_expr(e);
    return r;
}

// pinned antiderivatives, proofs of non-elementarity, unsupported input, and the parser
static void test_reg_integrate_exp(void) {
    const char* ok[][2] = {
        { "x*exp(x)", "(x - 1) * exp(x)" },
        { "x^2*exp(3x)", "(1/3x^2 - 2/9x + 2/27) * exp(3x)" },
        { "2x*exp(x^2)", "exp(x^2)" },
        { "(2x^2+1)*exp(x^2)", "x * exp(x^2)" },
        { "(x-1)*exp(x)/x^2", "exp(x) / x" },                     // a pole in the answer
        { "exp(-x)", "-exp(-x)" },
        { "exp(x)^2+exp(-x)", "-exp(-x) + 1/2 * exp(2x)" },
        { "(1+x)exp(x)+1/x", "x * exp(x) + ln(abs(x))" },          // the θ-free part uses LRT/Hermite
        { "exp(2x)*exp(x/2)", "2/5 * exp(5/2x)" },                 // one tower, θ = exp(x/2)
        { "exp(x)/exp(x)", "x" },
        { "x^3*exp(-x^2)", "(-1/2x^2 - 1/2) * exp(-x^2)" },
        { "-exp(-x)+2exp(x)", "exp(-x) + 2 * exp(x)" },            // unary minus before a function
    };
    for (size_t i = 0; i < sizeof ok / sizeof ok[0]; i++) {
        RischResult r = risch_str(ok[i][0]);
        assert(r.status == RISCH_OK);
        expect_str(expr_to_string(r.result), ok[i][1]);
        free_expr(r.result);
    }
    const char* none[] = { "exp(x^2)", "exp(x)/x", "exp(x)/x^2", "exp(x)/(x^2+1)", "(x+1)*exp(x^2)" };
    for (size_t i = 0; i < sizeof none / sizeof none[0]; i++) {
        RischResult r = risch_str(none[i]);
        assert(r.status == RISCH_NONELEMENTARY && r.result);
        free_expr(r.result);
    }
    const char* outside[] = { "1/(1+exp(x))", "exp(x)*exp(x^2)", "exp(1/x)", "ln(x)", "exp(3)", "exp(exp(x))", "exp(x)^(1/2)", "(1+exp(x))^(-1)" };
    for (size_t i = 0; i < sizeof outside / sizeof outside[0]; i++) {
        RischResult r = risch_str(outside[i]);
        assert(r.status == RISCH_UNSUPPORTED && !r.result && r.reason);
    }
    const char* parsed[][2] = {                                    // parse, then canonical form
        { "x exp(x)", "x * exp(x)" }, { "2exp(x)", "2 * exp(x)" }, { "exp (x)^2", "exp(2x)" },
        { "-exp(-x)", "-exp(-x)" }, { "exp(0)*x", "x" },
    };
    for (size_t i = 0; i < sizeof parsed / sizeof parsed[0]; i++) {
        Expr* e = parse_input(parsed[i][0]);
        const char* why = NULL;
        Expr* c = simplify_exp(e, &why);
        assert(c);
        expect_str(expr_to_string(c), parsed[i][1]);
        free_expr(e); free_expr(c);
    }
    printf("  test_reg_integrate_exp PASSED (%zu exact, %zu non-elementary, %zu unsupported)\n",
           sizeof ok / sizeof ok[0], sizeof none / sizeof none[0], sizeof outside / sizeof outside[0]);
}

// a random rational function with denominator (x - r)^m or 1
static Expr* random_coeff(void) {
    long nc[4];
    int nd = (int)rnd(0, 3);
    for (int i = 0; i <= nd; i++) nc[i] = rnd(-5, 5);
    Expr* num = create_poly_expr(poly_si(nc, nd));
    if (rnd(0, 2) == 0) return num;
    Polynomial* lin = poly_si((long[]){ -rnd(-3, 3), 1 }, 1);
    Polynomial* den = poly_pow(lin, (unsigned)rnd(1, 2));
    free_polynomial(lin);
    return create_div(num, create_poly_expr(den));
}

// F = Σ y_k·exp(k·v) with k in [-2, 2] \ {0}; v = c·x or c·x² + d·x (returned in *vout)
static Expr* random_exp_sum(Polynomial** vout) {
    long vc[3] = { 0, rnd(1, 3) * (rnd(0, 1) ? 1 : -1), rnd(0, 1) ? rnd(-2, 2) : 0 };
    Polynomial* v = poly_si(vc, vc[2] ? 2 : 1);
    *vout = poly_copy(v);
    Expr* F = NULL;
    for (int k = -2; k <= 2; k++) {
        if (k == 0 || rnd(0, 1)) continue;
        Rational kk = rat_from_int(k);
        Expr* t = create_mul(random_coeff(), create_exp(create_poly_expr(poly_scale(v, &kk))));
        rat_clear(&kk);
        F = F ? create_add(F, t) : t;
    }
    if (!F) F = create_mul(random_coeff(), create_exp(create_poly_expr(poly_copy(v))));
    free_polynomial(v);
    return F;
}

// equal as functions of x: the 5-point derivative of G against f, at points where two step
// sizes agree (a pole nearby makes the stencil unreliable, so such points are skipped). This
// check does not share the tower, so a wrong construction of v would show.
static bool derivative_matches(const Expr* G, const Expr* f) {
    const double xs[] = { -2.35, -1.85, -1.45, -1.15, -0.65, -0.35, 0.3, 0.6, 1.4, 1.7, 2.3, 0.85 };
    int checked = 0;                                       // poles sit at integers; exp(c·x²) grows fast
    for (size_t i = 0; i < sizeof xs / sizeof xs[0]; i++) {
        bool ok = true;
        double x0 = xs[i], v = expr_eval(f, x0, &ok), d1 = stencil(G, (Stencil){ x0, 1e-3 }, &ok),
               d2 = stencil(G, (Stencil){ x0, 5e-4 }, &ok);
        if (!ok || !isfinite(v) || fabs(v) > 1e8 || fabs(d1 - d2) > 1e-9 * (1 + fabs(v))) continue;
        if (fabs(d2 - v) > 1e-7 * (1 + fabs(v))) return false;
        checked++;
    }
    return checked >= 2;
}

// round trip: integrate(d/dx F) gives F's exp terms back exactly (the y_k are unique), and
// adding c·exp(k·v)/(x - r) to a solvable integrand makes it non-elementary, also when r is
// already a double pole of that coefficient (so the remainder test, not the simple-pole test, decides)
static void test_prop_integrate_exp(void) {
    int roundtrips = 0, rejected = 0;
    for (int trial = 0; trial < TRIALS; trial++) {
        Polynomial* v;
        Expr* F = random_exp_sum(&v);
        const char* why = NULL;
        Expr* f = differentiate_exp(F, &why);
        assert(f);
        RischResult r = integrate_exp(f);
        assert(r.status == RISCH_OK);
        const char* w1 = NULL;                                     // canonical forms spell out each
        Expr *cF = simplify_exp(F, &w1), *cI = simplify_exp(r.result, &w1); // exp(k·v), whatever the tower
        assert(cF && cI);
        expect_str(expr_to_string(cI), expr_to_string_free(cF));
        free_expr(cI);
        assert(derivative_matches(r.result, f));
        roundtrips++;
        long root = trial % 2 ? 1 : rnd(-3, 3);                   // odd trials: on a double pole of y = 1/(x-1)
        Expr* G = create_add(expr_copy(F), create_div(create_exp(create_poly_expr(poly_copy(v))), create_poly_expr(poly_si((long[]){ -1, 1 }, 1))));
        Expr* g = differentiate_exp(G, &why);
        assert(g);
        Expr* bad = create_add(expr_copy(g), create_div(create_mul(create_const(rnd(1, 5)), create_exp(create_poly_expr(poly_copy(v)))),
                                                         create_poly_expr(poly_si((long[]){ -root, 1 }, 1))));
        RischResult rb = integrate_exp(bad);
        assert(rb.status == RISCH_NONELEMENTARY);                   // always in the setting now
        rejected++;
        free_expr(rb.result);
        free_expr(G); free_expr(g); free_expr(bad);
        free_expr(r.result); free_expr(f); free_expr(F); free_polynomial(v);
    }
    printf("  test_prop_integrate_exp PASSED (%d round trips, %d non-elementary by construction)\n", roundtrips, rejected);
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
    test_reg_refine_shared_endpoint();
    test_reg_integrate_complete();
    test_reg_integrate_lrt();
    test_reg_integrate_radicals();
    test_reg_integrate_definite();
    test_reg_expr_lead();
    test_reg_qsurd();
    test_reg_radical_products();
    test_reg_integrate_exp();
    test_lrt_direct();
    printf("\nProperty Tests (seed %s):\n", env ? env : "default");
    test_prop_factor_roundtrip();
    test_prop_root_isolation();
    test_prop_integrate_derivative();
    test_prop_linear_systems();
    test_prop_print_parse_roundtrip();
    test_prop_subres_chain();
    test_prop_integrate_radicals();
    test_prop_exact_roots();
    test_prop_integrate_definite();
    test_prop_integrate_infinite();
    test_prop_definite_symmetric();
    test_prop_modular_gcd();
    test_prop_integrate_exp();
}
