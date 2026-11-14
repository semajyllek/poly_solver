// Subresultant chain of d(x) and g(t, x) = a(x) - t·d'(x) over Q[t], by the subresultant
// pseudo-remainder sequence with Lazard's shortcut for defective blocks (Ducos 2000).
//
// Invariant of the loop: A = S_d is a regular subresultant of degree d with principal
// coefficient s, and B = S_{d-1} has some degree e <= d - 1. Then, with δ = d - e:
//   S_e     = (lc(B) / s)^(δ-1) · B                 (zero for e < j < d - 1)
//   S_{e-1} = prem(A, -B) / (s^δ · lc(A))
// Every division is exact in Q[t]. The start is A = S_{m-1} = g, s = lc(g), B = prem(d, -g).
#include "subres.h"
#include <assert.h>
#include <stdlib.h>

// --- bivariate lifecycle ---

PolyTX ptx_alloc(int deg) {
    PolyTX a = { malloc((size_t)(deg + 1) * sizeof(Polynomial*)), deg };
    for (int k = 0; k <= deg; k++) a.c[k] = create_polynomial();
    return a;
}

PolyTX ptx_copy(const PolyTX* a) {
    PolyTX r = { malloc((size_t)(a->deg + 1) * sizeof(Polynomial*)), a->deg };
    for (int k = 0; k <= a->deg; k++) r.c[k] = poly_copy(a->c[k]);
    return r;
}

void ptx_free(PolyTX* a) {
    for (int k = 0; k <= a->deg; k++) free_polynomial(a->c[k]);
    free(a->c);
    a->c = NULL;
    a->deg = -1;
}

bool ptx_is_zero(const PolyTX* a) {
    for (int k = 0; k <= a->deg; k++)
        if (!poly_is_zero(a->c[k])) return false;
    return true;
}

// --- bivariate arithmetic (coefficients in Q[t]) ---

// replace coefficient k with an owned polynomial
static void set_coeff(PolyTX* a, int k, Polynomial* v) {
    free_polynomial(a->c[k]);
    a->c[k] = v;
}

// drop zero leading x-coefficients, keeping deg >= 0
static void ptx_trim(PolyTX* a) {
    while (a->deg > 0 && poly_is_zero(a->c[a->deg])) {
        free_polynomial(a->c[a->deg]);
        a->deg--;
    }
}

static const Polynomial* ptx_lc(const PolyTX* a) { return a->c[a->deg]; }

// a <- f·a
static void ptx_mul(PolyTX* a, const Polynomial* f) {
    for (int k = 0; k <= a->deg; k++) set_coeff(a, k, poly_multiply(a->c[k], f));
}

// a <- a / f, every coefficient division exact
static void ptx_div(PolyTX* a, const Polynomial* f) {
    for (int k = 0; k <= a->deg; k++) {
        PolyDivResult qr = poly_divmod(a->c[k], f);
        assert(poly_is_zero(qr.remainder)); // the chain identities guarantee exactness
        set_coeff(a, k, qr.quotient);
        free_polynomial(qr.remainder);
    }
}

static void ptx_neg(PolyTX* a) {
    Rational m1 = rat_from_int(-1);
    for (int k = 0; k <= a->deg; k++) set_coeff(a, k, poly_scale(a->c[k], &m1));
    rat_clear(&m1);
}

// r <- lc(b)·r - lc(r)·x^(deg r - deg b)·b, cancelling the leading term of r
static void prem_step(PolyTX* r, const PolyTX* b) {
    int shift = r->deg - b->deg;
    Polynomial* lr = poly_copy(ptx_lc(r));
    ptx_mul(r, ptx_lc(b));
    for (int k = 0; k <= b->deg; k++) {
        Polynomial* t = poly_multiply(lr, b->c[k]);
        set_coeff(r, k + shift, poly_subtract(r->c[k + shift], t));
        free_polynomial(t);
    }
    free_polynomial(lr);
    ptx_trim(r);
}

// prem(a, -b) = (-1)^(deg a - deg b + 1) · lc(b)^(deg a - deg b + 1) · a mod b
static PolyTX prem_neg(const PolyTX* a, const PolyTX* b) {
    int n = a->deg - b->deg + 1, steps = n;
    PolyTX r = ptx_copy(a);
    while (!ptx_is_zero(&r) && r.deg >= b->deg) { prem_step(&r, b); steps--; }
    Polynomial* f = poly_pow(ptx_lc(b), (unsigned)steps); // the multiplications skipped
    ptx_mul(&r, f);
    free_polynomial(f);
    if (n % 2) ptx_neg(&r);
    return r;
}

// --- the two ends of the chain ---

// d(x) as a bivariate polynomial constant in t
static PolyTX ptx_from_x(const Polynomial* d) {
    int m = polynomial_degree(d);
    PolyTX s = ptx_alloc(m);
    Rational c;
    rat_init(&c);
    for (int k = 0; k <= m; k++) { poly_coeff(&c, d, k); add_term(s.c[k], &c, 0); }
    rat_clear(&c);
    return s;
}

PolyTX rt_g(const RTPair* p) {
    int n = polynomial_degree(p->d) - 1;
    Polynomial* dp = poly_derivative(p->d);
    PolyTX g = ptx_alloc(n);
    Rational c;
    rat_init(&c);
    for (int k = 0; k <= n; k++) { // coefficient of x^k: a_k - d'_k · t
        poly_coeff(&c, p->a, k);
        add_term(g.c[k], &c, 0);
        poly_coeff(&c, dp, k);
        rat_neg(&c, &c);
        add_term(g.c[k], &c, 1);
    }
    rat_clear(&c);
    free_polynomial(dp);
    return g;
}

// --- the chain ---

typedef struct {
    PolyTX A;       // S_d, regular
    PolyTX B;       // S_{d-1}
    Polynomial* s;  // principal coefficient of S_d
} ChainStep;

// S_e = (lc(B) / s)^(δ-1) · B (Lazard); just B when the block is regular
static PolyTX lazard(const ChainStep* st) {
    unsigned k = (unsigned)(st->A.deg - st->B.deg - 1);
    PolyTX r = ptx_copy(&st->B);
    if (k == 0) return r;
    Polynomial* num = poly_pow(ptx_lc(&st->B), k);
    Polynomial* den = poly_pow(st->s, k);
    ptx_mul(&r, num);
    ptx_div(&r, den);
    free_polynomial(num); free_polynomial(den);
    return r;
}

// S_{e-1} = prem(A, -B) / (s^δ · lc(A))
static PolyTX next_sub(const ChainStep* st) {
    PolyTX r = prem_neg(&st->A, &st->B);
    Polynomial* sd = poly_pow(st->s, (unsigned)(st->A.deg - st->B.deg));
    Polynomial* den = poly_multiply(sd, ptx_lc(&st->A));
    ptx_div(&r, den);
    free_polynomial(sd); free_polynomial(den);
    return r;
}

static void put(SubresChain* ch, int j, PolyTX v) {
    ptx_free(&ch->S[j]);
    ch->S[j] = v;
}

// record S_{d-1} and S_e, then move to (S_e, S_{e-1}); false once the chain is complete
static bool advance(SubresChain* ch, ChainStep* st) {
    int d = st->A.deg, e = st->B.deg;
    put(ch, d - 1, ptx_copy(&st->B));
    PolyTX C = lazard(st);
    if (e < d - 1) put(ch, e, ptx_copy(&C));
    if (e == 0) { ptx_free(&C); return false; }
    PolyTX next = next_sub(st);
    ptx_free(&st->A); ptx_free(&st->B); free_polynomial(st->s);
    st->s = poly_copy(ptx_lc(&C));
    st->A = C;
    st->B = next;
    return true;
}

SubresChain rt_chain(const RTPair* p) {
    int m = polynomial_degree(p->d);
    SubresChain ch = { malloc((size_t)(m + 1) * sizeof(PolyTX)), m };
    for (int j = 0; j <= m; j++) ch.S[j] = ptx_alloc(0);
    put(&ch, m, ptx_from_x(p->d));
    put(&ch, m - 1, rt_g(p));
    ChainStep st = { ptx_copy(&ch.S[m - 1]), prem_neg(&ch.S[m], &ch.S[m - 1]),
                     poly_copy(ptx_lc(&ch.S[m - 1])) };
    while (!ptx_is_zero(&st.B) && advance(&ch, &st)) {}
    ptx_free(&st.A); ptx_free(&st.B); free_polynomial(st.s);
    return ch;
}

void chain_free(SubresChain* ch) {
    for (int j = 0; j <= ch->m; j++) ptx_free(&ch->S[j]);
    free(ch->S);
    ch->S = NULL;
}
