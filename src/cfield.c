// Q(i)[x] as pairs of rational polynomials, Q(i)(x) with a common real denominator, and the
// Risch differential equation y' + f·y = g over Q(i)(x) (the algorithm of risch.c's solve_rde,
// with complex leading coefficients).
#include "cfield.h"
#include "modgcd.h"

static Polynomial* poly_one(void) {
    Rational one = rat_one_val();
    Polynomial* p = poly_constant(&one);
    rat_clear(&one);
    return p;
}

// --- Q(i)[x] ---

CPoly cp_zero(void) { return (CPoly){ create_polynomial(), create_polynomial() }; }
CPoly cp_copy(const CPoly* a) { return (CPoly){ poly_copy(a->re), poly_copy(a->im) }; }
void cp_free(CPoly* a) { free_polynomial(a->re); free_polynomial(a->im); a->re = a->im = NULL; }
bool cp_is_zero(const CPoly* a) { return poly_is_zero(a->re) && poly_is_zero(a->im); }

int cp_deg(const CPoly* a) {
    int r = poly_is_zero(a->re) ? -1 : polynomial_degree(a->re), i = poly_is_zero(a->im) ? -1 : polynomial_degree(a->im);
    return r > i ? r : i;
}

CPoly cp_add(const CPoly* a, const CPoly* b, int sign) {
    if (sign > 0) return (CPoly){ poly_add(a->re, b->re), poly_add(a->im, b->im) };
    return (CPoly){ poly_subtract(a->re, b->re), poly_subtract(a->im, b->im) };
}

bool cp_equal(const CPoly* a, const CPoly* b) {
    CPoly d = cp_add(a, b, -1);
    bool z = cp_is_zero(&d);
    cp_free(&d);
    return z;
}

// (a + bi)(c + di) = (ac - bd) + (ad + bc)i
CPoly cp_mul(const CPoly* a, const CPoly* b) {
    Polynomial *ac = poly_multiply(a->re, b->re), *bd = poly_multiply(a->im, b->im);
    Polynomial *ad = poly_multiply(a->re, b->im), *bc = poly_multiply(a->im, b->re);
    CPoly r = { poly_subtract(ac, bd), poly_add(ad, bc) };
    free_polynomial(ac); free_polynomial(bd); free_polynomial(ad); free_polynomial(bc);
    return r;
}

CPoly cp_neg(const CPoly* a) {
    CPoly z = cp_zero(), r = cp_add(&z, a, -1);
    cp_free(&z);
    return r;
}

CPoly cp_derivative(const CPoly* a) { return (CPoly){ poly_derivative(a->re), poly_derivative(a->im) }; }

// a·p for a real polynomial p
static CPoly cp_times_real(const CPoly* a, const Polynomial* p) {
    return (CPoly){ poly_multiply(a->re, p), poly_multiply(a->im, p) };
}

// --- Q(i)(x) ---

CFn cfn_from(CPoly n, Polynomial* den) {
    if (cp_is_zero(&n)) { cp_free(&n); free_polynomial(den); return (CFn){ cp_zero(), poly_one() }; }
    Polynomial* g0 = poly_gcd_fast(n.re, n.im);
    Polynomial* g = poly_gcd_fast(g0, den);               // monic; gcd(0, den) = den/lc
    Polynomial *re = poly_divide(n.re, g), *im = poly_divide(n.im, g), *d = poly_divide(den, g);
    Rational lc;                                          // monic denominator
    rat_init(&lc);
    poly_coeff(&lc, d, polynomial_degree(d));
    mpq_inv(lc.val, lc.val);
    CFn r = { { poly_scale(re, &lc), poly_scale(im, &lc) }, poly_scale(d, &lc) };
    free_polynomial(re); free_polynomial(im); free_polynomial(d); free_polynomial(g0); free_polynomial(g);
    cp_free(&n); free_polynomial(den);
    rat_clear(&lc);
    return r;
}

CFn cfn_real(const Polynomial* num, const Polynomial* den) {
    return cfn_from((CPoly){ poly_copy(num), create_polynomial() }, poly_copy(den));
}

CFn cfn_const(long re_num, long im_num, long den) {
    Rational a, b;
    rat_init(&a); rat_init(&b);
    rat_set_si(&a, re_num, den);
    rat_set_si(&b, im_num, den);
    CFn r = cfn_from((CPoly){ poly_constant(&a), poly_constant(&b) }, poly_one());
    rat_clear(&a); rat_clear(&b);
    return r;
}

CFn cfn_poly(const CPoly* a) { return cfn_from(cp_copy(a), poly_one()); }

CFn cfn_zero(void) { return (CFn){ cp_zero(), poly_one() }; }
CFn cfn_copy(const CFn* a) { return (CFn){ cp_copy(&a->n), poly_copy(a->den) }; }
void cfn_free(CFn* a) { cp_free(&a->n); free_polynomial(a->den); a->den = NULL; }
bool cfn_is_zero(const CFn* a) { return cp_is_zero(&a->n); }

CFn cfn_add(const CFn* a, const CFn* b, int sign) {
    CPoly l = cp_times_real(&a->n, b->den), r = cp_times_real(&b->n, a->den);
    CPoly s = cp_add(&l, &r, sign);
    cp_free(&l); cp_free(&r);
    return cfn_from(s, poly_multiply(a->den, b->den));
}

CFn cfn_mul(const CFn* a, const CFn* b) { return cfn_from(cp_mul(&a->n, &b->n), poly_multiply(a->den, b->den)); }

// d/(re + i·im) = d·(re - i·im)/(re² + im²)
CFn cfn_inv(const CFn* a) {
    Polynomial *r2 = poly_multiply(a->n.re, a->n.re), *i2 = poly_multiply(a->n.im, a->n.im);
    CPoly conj = { poly_copy(a->n.re), poly_copy(a->n.im) };
    Rational m1 = rat_from_int(-1);
    Polynomial* nim = poly_scale(conj.im, &m1);
    free_polynomial(conj.im);
    conj.im = nim;
    CPoly n = cp_times_real(&conj, a->den);
    CFn r = cfn_from(n, poly_add(r2, i2));
    cp_free(&conj); free_polynomial(r2); free_polynomial(i2); rat_clear(&m1);
    return r;
}

// (n/d)' = (n'·d - n·d')/d²
CFn cfn_derivative(const CFn* a) {
    CPoly dn = cp_derivative(&a->n), l = cp_times_real(&dn, a->den);
    Polynomial* dd = poly_derivative(a->den);
    CPoly r = cp_times_real(&a->n, dd), s = cp_add(&l, &r, -1);
    cp_free(&dn); cp_free(&l); cp_free(&r); free_polynomial(dd);
    return cfn_from(s, poly_multiply(a->den, a->den));
}

// --- the Risch differential equation ---

// a + bi as two rationals
typedef struct {
    Rational re, im;
} CRat;

static CRat cp_coeff(const CPoly* a, int k) {
    CRat c;
    rat_init(&c.re); rat_init(&c.im);
    poly_coeff(&c.re, a->re, k);
    poly_coeff(&c.im, a->im, k);
    return c;
}

static void crat_clear(CRat* c) { rat_clear(&c->re); rat_clear(&c->im); }

// (a + bi)/(c + di) = ((ac + bd) + (bc - ad)i)/(c² + d²)
static CRat crat_div(const CRat* x, const CRat* y) {
    CRat q;
    Rational t, n;
    rat_init(&q.re); rat_init(&q.im); rat_init(&t); rat_init(&n);
    rat_mul(&n, &y->re, &y->re);
    rat_mul(&t, &y->im, &y->im);
    rat_add(&n, &n, &t);
    rat_mul(&q.re, &x->re, &y->re);
    rat_mul(&t, &x->im, &y->im);
    rat_add(&q.re, &q.re, &t);
    rat_mul(&q.im, &x->im, &y->re);
    rat_mul(&t, &x->re, &y->im);
    rat_sub(&q.im, &q.im, &t);
    rat_div(&q.re, &q.re, &n);
    rat_div(&q.im, &q.im, &n);
    rat_clear(&t); rat_clear(&n);
    return q;
}

// c·x^k
static CPoly cp_monomial(const CRat* c, int k) {
    CPoly t = cp_zero();
    add_term(t.re, &c->re, k);
    add_term(t.im, &c->im, k);
    return t;
}

// q with A·q' + B·q = C (A real), deg B >= deg A, triangular from the top; false if none
static bool c_solve_poly_rde(const Polynomial* A, const CPoly* B, CPoly* C) {
    int db = cp_deg(B);
    CRat lb = cp_coeff(B, db);
    CPoly q = cp_zero();
    while (!cp_is_zero(C) && cp_deg(C) >= db) {
        int dc = cp_deg(C);
        CRat lr = cp_coeff(C, dc), c = crat_div(&lr, &lb);
        CPoly t = cp_monomial(&c, dc - db), dt = cp_derivative(&t);
        CPoly At = cp_times_real(&dt, A), Bt = cp_mul(B, &t), s = cp_add(&At, &Bt, 1);
        CPoly nr = cp_add(C, &s, -1), nq = cp_add(&q, &t, 1);
        cp_free(C); cp_free(&q);
        *C = nr; q = nq;
        cp_free(&t); cp_free(&dt); cp_free(&At); cp_free(&Bt); cp_free(&s);
        crat_clear(&lr); crat_clear(&c);
    }
    crat_clear(&lb);
    bool solved = cp_is_zero(C);
    cp_free(C);
    *C = q;                                               // the solution, when solved
    return solved;
}

// poles of y only at poles of g, one order higher (f has none), so no simple poles; then
// y = q/h with h = gcd(D, D') and h·q' + (f·h - h')·q = N·h²/D
bool c_solve_rde(const CPoly* f, const CFn* g, CFn* y) {
    if (cfn_is_zero(g)) { *y = cfn_zero(); return true; }
    Polynomial *Dp = poly_derivative(g->den), *h = poly_gcd_fast(g->den, Dp), *h2 = poly_multiply(h, h);
    PolyDivResult qd = poly_divmod(h2, g->den);           // D | h² iff every pole has order >= 2
    bool ok = poly_is_zero(qd.remainder);
    if (ok) {
        CPoly C = cp_times_real(&g->n, qd.quotient), fh = cp_times_real(f, h);
        Polynomial* hp = poly_derivative(h);
        Polynomial* bre = poly_subtract(fh.re, hp);
        free_polynomial(fh.re);
        fh.re = bre;                                      // B = f·h - h'
        ok = c_solve_poly_rde(h, &fh, &C);
        if (ok) *y = cfn_from(C, poly_copy(h));
        else cp_free(&C);
        cp_free(&fh); free_polynomial(hp);
    }
    free_polynomial(Dp); free_polynomial(h); free_polynomial(h2);
    free_polynomial(qd.quotient); free_polynomial(qd.remainder);
    return ok;
}
