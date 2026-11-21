// Test oracle: Collins subresultants of d and g = a - t·d' straight from the determinant
// definition. With F = d of degree m and G = g of formal degree n = m - 1, the matrix M_j
// stacks the n-j shifted copies x^k·F and the m-j shifted copies x^k·G; the coefficient of
// x^i in S_j is the determinant of the first m+n-2j-1 columns of M_j plus the column of x^i.
// Every entry is linear in t, so each determinant has degree <= m - j in t: we evaluate it
// exactly at m - j + 1 integer points and interpolate. Slow, but independent of the chain.
#include "test_subres_ref.h"
#include "matrix.h"
#include <stdlib.h>

// --- dense coefficient arrays (index = power of x) ---

static Rational* dense(const Polynomial* p, int n) {
    Rational* c = malloc((size_t)(n + 1) * sizeof(Rational));
    for (int k = 0; k <= n; k++) { rat_init(&c[k]); poly_coeff(&c[k], p, k); }
    return c;
}

static void free_dense(Rational* c, int n) {
    for (int k = 0; k <= n; k++) rat_clear(&c[k]);
    free(c);
}

// g(t0, x) = a(x) - t0·d'(x) with formal degree m - 1
static Rational* eval_g(const RTPair* p, const Rational* t0) {
    int n = polynomial_degree(p->d) - 1;
    Polynomial* dp = poly_derivative(p->d);
    Rational* a = dense(p->a, n);
    Rational* d = dense(dp, n);
    Rational tmp;
    rat_init(&tmp);
    for (int k = 0; k <= n; k++) {
        rat_mul(&tmp, t0, &d[k]);
        rat_sub(&a[k], &a[k], &tmp);
    }
    rat_clear(&tmp);
    free_dense(d, n);
    free_polynomial(dp);
    return a;
}

// --- Sylvester submatrix determinants ---

typedef struct {
    const Rational* f; int m; // F coefficients, degree m
    const Rational* g; int n; // G coefficients, formal degree n
} SylvPair;

typedef struct {
    int j; // which subresultant
    int i; // which coefficient (power of x)
} SubresIdx;

// one row of the Sylvester matrix: x^shift times a dense polynomial of degree deg
typedef struct {
    const Rational* c;
    int deg;
    int shift;
} SylvRow;

// coefficient of x^power in the row's polynomial
static const Rational* row_coeff(const SylvRow* row, int power) {
    int k = power - row->shift;
    return (k >= 0 && k <= row->deg) ? &row->c[k] : NULL;
}

// coefficient of x^i in S_j, as an exact determinant
static void subres_det(Rational* out, const SylvPair* s, const SubresIdx* idx) {
    int j = idx->j, L = s->m + s->n - 2 * j, nf = s->n - j;
    Matrix* M = matrix_create(L, L);
    for (int r = 0; r < L; r++) {
        bool is_f = r < nf;
        SylvRow row = is_f ? (SylvRow){ s->f, s->m, nf - 1 - r }
                           : (SylvRow){ s->g, s->n, (s->m - j - 1) - (r - nf) };
        for (int c = 0; c < L; c++) {
            int power = (c < L - 1) ? s->m + s->n - j - 1 - c : idx->i;
            const Rational* v = row_coeff(&row, power);
            if (v) rat_set(&MAT_AT(M, r, c), v);
        }
    }
    matrix_det(out, M);
    matrix_free(M);
}

// --- interpolation ---

static Polynomial* interpolate(const Rational* ts, const Rational* vs, int n) {
    Rational* c = malloc((size_t)n * sizeof(Rational));
    Rational num, den;
    rat_init(&num); rat_init(&den);
    for (int i = 0; i < n; i++) rat_init_set(&c[i], &vs[i]);
    for (int k = 1; k < n; k++) // Newton divided differences, in place
        for (int i = n - 1; i >= k; i--) {
            rat_sub(&num, &c[i], &c[i - 1]);
            rat_sub(&den, &ts[i], &ts[i - k]);
            rat_div(&c[i], &num, &den);
        }
    Polynomial* p = poly_constant(&c[n - 1]);
    Rational one = rat_one_val();
    for (int k = n - 2; k >= 0; k--) { // Horner: p = p·(t - ts[k]) + c[k]
        Polynomial* lin = create_polynomial();
        add_term(lin, &one, 1);
        rat_neg(&num, &ts[k]);
        add_term(lin, &num, 0);
        Polynomial* q = poly_multiply(p, lin);
        add_term(q, &c[k], 0);
        free_polynomial(p); free_polynomial(lin);
        p = q;
    }
    for (int i = 0; i < n; i++) rat_clear(&c[i]);
    free(c);
    rat_clear(&num); rat_clear(&den); rat_clear(&one);
    return p;
}

// --- resultant and subresultants in t ---

// values[i][v] = coefficient of x^i of S_j at t = v, for v = 0..npts-1
static void sample_subres(const RTPair* p, int j, Rational** values) {
    int m = polynomial_degree(p->d), npts = m - j + 1;
    Rational* f = dense(p->d, m);
    for (int v = 0; v < npts; v++) {
        Rational t0 = rat_from_int(v);
        Rational* g = eval_g(p, &t0);
        SylvPair s = { f, m, g, m - 1 };
        for (int i = 0; i <= j; i++) {
            SubresIdx idx = { j, i };
            subres_det(&values[i][v], &s, &idx);
        }
        free_dense(g, m - 1);
        rat_clear(&t0);
    }
    free_dense(f, m);
}

// S_j(t, x) for 0 <= j <= deg d - 2, its coefficients as polynomials in t
PolyTX ref_subresultant(const RTPair* p, int j) {
    int m = polynomial_degree(p->d), npts = m - j + 1;
    Rational** values = malloc((size_t)(j + 1) * sizeof(Rational*));
    Rational* ts = malloc((size_t)npts * sizeof(Rational));
    for (int v = 0; v < npts; v++) { rat_init(&ts[v]); mpq_set_si(ts[v].val, v, 1); }
    for (int i = 0; i <= j; i++) {
        values[i] = malloc((size_t)npts * sizeof(Rational));
        for (int v = 0; v < npts; v++) rat_init(&values[i][v]);
    }
    sample_subres(p, j, values);
    PolyTX S = { malloc((size_t)(j + 1) * sizeof(Polynomial*)), j };
    for (int i = 0; i <= j; i++) {
        S.c[i] = interpolate(ts, values[i], npts);
        for (int v = 0; v < npts; v++) rat_clear(&values[i][v]);
        free(values[i]);
    }
    for (int v = 0; v < npts; v++) rat_clear(&ts[v]);
    free(ts);
    free(values);
    return S;
}

