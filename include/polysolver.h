// polysolver: exact symbolic algebra over the rationals (GMP). The whole public API, in one header.
//
//   cc my_prog.c $(pkg-config --cflags --libs polysolver)
//
// Every name is prefixed ps_ (PS_ for macros and enum constants); the library defines nothing else
// that a program can collide with. Ownership: functions return new allocations that the caller frees
// (ps_expr_free, ps_poly_free, free() for char* strings, ...). The ps_expr_ constructors take
// ownership of their ps_expr and ps_poly arguments; every other function leaves its inputs alone.
#ifndef POLYSOLVER_H
#define POLYSOLVER_H

#include <gmp.h>
#include <stdbool.h>

#define PS_VERSION_MAJOR 0
#define PS_VERSION_MINOR 1
#define PS_VERSION_PATCH 0
#define PS_VERSION "0.1.0"

#if defined(__GNUC__) || defined(__clang__)
#define PS_API __attribute__((visibility("default")))
#else
#define PS_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

// --- types ---

typedef struct {                   // an exact rational number (an mpq_t; init before use, clear after)
    mpq_t val;
} ps_rational;

typedef struct ps_poly ps_poly;    // a polynomial in x with rational coefficients (opaque)
typedef struct ps_expr ps_expr;    // an expression tree in x (opaque)

typedef struct {                   // ps_poly_divmod: dividend = quotient·divisor + remainder
    ps_poly* quotient;
    ps_poly* remainder;
} ps_divmod_result;

typedef struct {                   // ps_poly_xgcd: s·a + t·b = g, g = gcd(a, b) monic
    ps_poly* s;
    ps_poly* t;
    ps_poly* g;
} ps_bezout;

typedef struct {                   // poly = content · Π factors[i]^multiplicities[i], factors primitive
    ps_poly** factors;
    int* multiplicities;
    int count;
    ps_rational content;
} ps_factorization;

typedef struct {                   // poly = Π parts[i]^(i+1), each part square-free
    ps_poly** parts;
    int count;
} ps_sqfree;

typedef struct {                   // [lo, hi] holding exactly one real root (lo = hi: the root itself)
    ps_rational lo, hi;
} ps_root_interval;

typedef struct {
    ps_root_interval* items;       // sorted ascending
    int count;
} ps_root_list;

typedef struct {
    ps_poly* factor;               // irreducible over Q
    int multiplicity;
    ps_root_interval* real;        // certified isolating intervals of its real roots, ascending
    int num_real;
} ps_irrational_factor;

typedef struct {
    ps_rational* rational_roots;   // exact, repeated by multiplicity
    int num_rational;
    ps_irrational_factor* irr;     // one entry per distinct irreducible factor of degree >= 2
    int num_irr;
    double* irrational_roots;      // approximations of the real irrational roots, by multiplicity
    int num_irrational;
    ps_poly** irreducible_factors; // factors that have non-real roots, by multiplicity
    int num_irreducible;
} ps_solutions;

typedef enum { PS_RISCH_OK, PS_RISCH_NONELEMENTARY, PS_RISCH_UNSUPPORTED } ps_risch_status;

typedef struct {
    ps_risch_status status;
    ps_expr* result;               // OK: the antiderivative; NONELEMENTARY: the part whose integral is not elementary
    const char* reason;            // otherwise: why (a static string)
} ps_risch_result;

typedef struct {
    ps_expr* exact;                // the exact value, or NULL when only the number is available
    double approx;                 // tanh-sinh or adaptive quadrature, independent of the antiderivative
    const char* error;             // set when there is no value at all (a pole, a divergent end)
    const char* note;              // set when exact is NULL but approx is valid: the reason
} ps_definite_result;

typedef struct {                   // a bound: the rational *x (x·π when pi), or ±∞ when inf = ±1
    const ps_rational* x;
    int inf;
    bool pi;
} ps_endpoint;

typedef struct {
    ps_endpoint a, b;
} ps_bounds;

typedef struct {                   // row-major: element (r, c) is PS_MAT_AT(m, r, c)
    ps_rational* data;
    int rows;
    int cols;
} ps_matrix;

#define PS_MAT_AT(m, r, c) ((m)->data[(r) * (m)->cols + (c)])

typedef struct {
    ps_rational* data;
    int n;
} ps_vector;

typedef enum { PS_LINSYS_UNIQUE, PS_LINSYS_INCONSISTENT, PS_LINSYS_INFINITE } ps_linsys_status;

typedef struct {
    ps_linsys_status status;
    ps_vector* solution;           // NULL if inconsistent
    int rank;
} ps_linsys_result;

// --- the library ---

PS_API const char* ps_version(void);                                               // PS_VERSION of the built library
PS_API const ps_poly* ps_expr_poly(const ps_expr* e);                              // e as a polynomial (borrowed), or NULL

// --- rationals ---

PS_API void ps_rat_init(ps_rational* r);
PS_API void ps_rat_clear(ps_rational* r);
PS_API void ps_rat_init_set(ps_rational* dst, const ps_rational* src);
PS_API void ps_rat_set(ps_rational* dst, const ps_rational* src);
PS_API void ps_rat_set_si(ps_rational* r, long num, long den);
PS_API ps_rational ps_rat_from_si(long num, long den);
PS_API ps_rational ps_rat_from_int(long n);
PS_API ps_rational ps_rat_zero(void);
PS_API ps_rational ps_rat_one(void);
PS_API void ps_rat_add(ps_rational* result, const ps_rational* a, const ps_rational* b);
PS_API void ps_rat_sub(ps_rational* result, const ps_rational* a, const ps_rational* b);
PS_API void ps_rat_mul(ps_rational* result, const ps_rational* a, const ps_rational* b);
PS_API void ps_rat_div(ps_rational* result, const ps_rational* a, const ps_rational* b);
PS_API void ps_rat_neg(ps_rational* result, const ps_rational* a);
PS_API void ps_rat_abs(ps_rational* result, const ps_rational* a);
PS_API void ps_rat_pow(ps_rational* result, const ps_rational* base, int exp);
PS_API bool ps_rat_eq(const ps_rational* a, const ps_rational* b);
PS_API int ps_rat_cmp(const ps_rational* a, const ps_rational* b);
PS_API bool ps_rat_is_zero(const ps_rational* a);
PS_API bool ps_rat_is_one(const ps_rational* a);
PS_API bool ps_rat_is_negative(const ps_rational* a);
PS_API double ps_rat_to_double(const ps_rational* r);
PS_API char* ps_rat_to_string(const ps_rational* r);                                 // caller frees
PS_API char* ps_rat_to_decimal(const ps_rational* r, int digits);
PS_API const char* ps_rat_parse(ps_rational* out, const char* s);

// --- polynomials ---

PS_API ps_poly* ps_poly_new(void);
PS_API void ps_poly_free(ps_poly* poly);
PS_API ps_poly* ps_poly_copy(const ps_poly* poly);
PS_API void ps_poly_add_term(ps_poly* poly, const ps_rational* coeff, int exponent);
PS_API ps_poly* ps_poly_constant(const ps_rational* c);                              // the constant polynomial c
PS_API int ps_poly_degree(const ps_poly* poly);
PS_API void ps_poly_coeff(ps_rational* out, const ps_poly* poly, int degree);        // 0 when absent
PS_API bool ps_poly_is_zero(const ps_poly* poly);
PS_API ps_poly* ps_poly_add(const ps_poly* poly1, const ps_poly* poly2);
PS_API ps_poly* ps_poly_sub(const ps_poly* poly1, const ps_poly* poly2);
PS_API ps_poly* ps_poly_mul(const ps_poly* poly1, const ps_poly* poly2);
PS_API ps_poly* ps_poly_scale(const ps_poly* p, const ps_rational* c);
PS_API ps_poly* ps_poly_pow(const ps_poly* base, unsigned exp);                      // binary exponentiation
PS_API ps_divmod_result ps_poly_divmod(const ps_poly* dividend, const ps_poly* divisor);
PS_API ps_poly* ps_poly_quo(const ps_poly* dividend, const ps_poly* divisor);
PS_API ps_poly* ps_poly_rem(const ps_poly* dividend, const ps_poly* divisor);
PS_API ps_poly* ps_poly_gcd(const ps_poly* poly1, const ps_poly* poly2);
PS_API ps_poly* ps_poly_lcm(const ps_poly* poly1, const ps_poly* poly2);
PS_API ps_bezout ps_poly_xgcd(const ps_poly* a, const ps_poly* b);
PS_API void ps_bezout_free(ps_bezout* bz);
PS_API void ps_poly_make_monic(ps_poly* poly);
PS_API ps_poly* ps_poly_derivative(const ps_poly* poly);
PS_API ps_poly* ps_poly_integral(const ps_poly* poly);
PS_API void ps_poly_eval(ps_rational* out, const ps_poly* poly, const ps_rational* x); // exact Horner
PS_API double ps_poly_eval_double(const ps_poly* poly, double x);
PS_API char* ps_poly_to_string(const ps_poly* poly);                                 // caller frees

// --- factoring, roots, solving ---

PS_API ps_factorization* ps_factorize(const ps_poly* poly);
PS_API void ps_factorization_free(ps_factorization* f);
PS_API char* ps_factorization_to_string(const ps_factorization* f);
PS_API ps_sqfree ps_square_free(const ps_poly* poly);
PS_API void ps_square_free_free(ps_sqfree* d);
PS_API ps_root_list ps_real_roots(const ps_poly* f);
PS_API void ps_root_refine(ps_root_interval* r, const ps_poly* f, const ps_rational* width);
PS_API void ps_root_list_free(ps_root_list* L);
PS_API ps_solutions* ps_solve(const ps_poly* poly);
PS_API char* ps_solutions_to_string(const ps_solutions* sol);
PS_API void ps_solutions_free(ps_solutions* sol);

// --- expressions ---

PS_API ps_expr* ps_expr_int(long value);                                             // exact integer constant
PS_API ps_expr* ps_expr_rational(const ps_rational* value);
PS_API ps_expr* ps_expr_decimal(const char* literal);                                // exact decimal literal ("1.25", "2e30")
PS_API ps_expr* ps_expr_var(const char* name);
PS_API ps_expr* ps_expr_from_poly(ps_poly* poly);
PS_API ps_expr* ps_expr_ratfn(ps_poly* num, ps_poly* den);
PS_API ps_expr* ps_expr_add(ps_expr* left, ps_expr* right);
PS_API ps_expr* ps_expr_sub(ps_expr* left, ps_expr* right);
PS_API ps_expr* ps_expr_mul(ps_expr* left, ps_expr* right);
PS_API ps_expr* ps_expr_div(ps_expr* numerator, ps_expr* denominator);
PS_API ps_expr* ps_expr_pow(ps_expr* base, ps_expr* exponent);
PS_API ps_expr* ps_expr_neg(ps_expr* operand);
PS_API ps_expr* ps_expr_ln(ps_expr* arg);
PS_API ps_expr* ps_expr_exp(ps_expr* arg);
PS_API ps_expr* ps_expr_sin(ps_expr* arg);
PS_API ps_expr* ps_expr_cos(ps_expr* arg);
PS_API ps_expr* ps_expr_atan(ps_expr* arg);
PS_API ps_expr* ps_expr_sqrt(ps_expr* arg);
PS_API ps_expr* ps_expr_pi(void);                                                    // π
PS_API ps_expr* ps_expr_copy(const ps_expr* e);                                      // deep copy
PS_API void ps_expr_free(ps_expr* expr);
PS_API char* ps_expr_to_string(const ps_expr* expr);
PS_API double ps_expr_eval(const ps_expr* expr, double x, bool* ok);
PS_API ps_expr* ps_canonicalize(ps_expr* expr);
PS_API ps_expr* ps_simplify(ps_expr* expr);
PS_API ps_expr* ps_factor(ps_expr* expr);
PS_API ps_expr* ps_differentiate(ps_expr* expr, const char* variable);
PS_API ps_expr* ps_integrate(ps_expr* expr, const char* variable);                   // integrate.c
PS_API ps_expr* ps_integrate_rootsum(ps_expr* expr, const char* variable);           // same, but RootSum is never expanded
PS_API const char* ps_integrate_error(void);                                         // reason for the last NULL from integrate
PS_API ps_expr* ps_parse(const char* input);                                         // the parse tree, unsimplified (NULL with a message on error)
PS_API ps_expr* ps_parse_simplify(const char* input);

// --- exp/ln/sin/cos, square roots, definite integrals ---

PS_API ps_risch_result ps_integrate_exp(const ps_expr* f);
PS_API ps_expr* ps_differentiate_exp(const ps_expr* f, const char** reason);
PS_API ps_expr* ps_simplify_exp(const ps_expr* f, const char** reason);
PS_API ps_risch_result ps_integrate_sqrt(const ps_expr* f);                          // RISCH_OK or RISCH_UNSUPPORTED (never a proof)
PS_API ps_expr* ps_differentiate_sqrt(const ps_expr* f, const char** why);           // NULL with *why on failure
PS_API ps_expr* ps_simplify_sqrt(const ps_expr* f, const char** why);
PS_API ps_definite_result ps_integrate_definite(ps_expr* f, const ps_bounds* ab);
PS_API ps_definite_result ps_integrate_definite_sqrt(const ps_expr* f, const ps_bounds* ab);
PS_API bool ps_parse_endpoint(ps_rational* out, ps_endpoint* end, const char* text);
PS_API ps_expr* ps_parse_integrand(const char* text);
PS_API void ps_definite_free(ps_definite_result* r);

// --- linear algebra ---

PS_API ps_matrix* ps_matrix_new(int rows, int cols);
PS_API ps_matrix* ps_matrix_copy(const ps_matrix* m);
PS_API void ps_matrix_free(ps_matrix* m);
PS_API ps_matrix* ps_matrix_parse(const char* str);
PS_API char* ps_matrix_to_string(const ps_matrix* m);
PS_API ps_matrix* ps_matrix_rref(const ps_matrix* m);
PS_API int ps_matrix_rank(const ps_matrix* m);
PS_API void ps_matrix_det(ps_rational* result, const ps_matrix* m);
PS_API ps_matrix* ps_matrix_inverse(const ps_matrix* m);
PS_API ps_vector* ps_vector_new(int n);
PS_API ps_vector* ps_vector_copy(const ps_vector* v);
PS_API void ps_vector_free(ps_vector* v);
PS_API ps_vector* ps_vector_parse(const char* str);
PS_API char* ps_vector_to_string(const ps_vector* v);
PS_API ps_linsys_result* ps_linsys_solve(const ps_matrix* A, const ps_vector* b);
PS_API char* ps_linsys_to_string(const ps_linsys_result* r);
PS_API void ps_linsys_free(ps_linsys_result* r);

#ifdef __cplusplus
}
#endif

#endif
