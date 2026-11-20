// K[θ] for K = Q(x): division with remainder, gcd, Bezout and square-free decomposition,
// all over exact QFn coefficients.
#include "kpoly.h"
#include "factor.h"
#include <stdlib.h>

static QFn qfn_one(void) {
    Rational one = rat_one_val();
    QFn r = { poly_constant(&one), poly_constant(&one) };
    rat_clear(&one);
    return r;
}

int kp_deg(const LPoly* a) { return a->n == 0 ? -1 : a->lo + a->n - 1; }

LPoly kp_one(void) { return lpoly_term(qfn_one(), 0); }

bool kp_is_one(const LPoly* a) {
    if (a->n != 1 || a->lo != 0) return false;
    QFn one = qfn_one(), d = qfn_add(&a->a[0], &one, -1);
    bool z = qfn_is_zero(&d);
    qfn_free(&one); qfn_free(&d);
    return z;
}

// --- the algorithms (kpoly_tmpl.h) over QFn ---

// k·a
static QFn qfn_scale_int(const QFn* a, int k) {
    Rational kk = rat_from_int(k);
    QFn c = { poly_scale(a->num, &kk), poly_copy(a->den) };
    rat_clear(&kk);
    return c;
}

#define P LPoly
#define C QFn
#define CTX int
#define PFX(name) kq_##name
#define KDIV KDiv
#define KBEZ KBezout
#define KSF KSquareFree
#define P_DEG(p) kp_deg(p)
#define P_LC(p) (&(p)->a[(p)->n - 1])
#define P_TERM(x, c, k) lpoly_term(qfn_copy(c), (k))
#define P_ADD(x, a, b, sign) lpoly_add((a), (b), (sign))
#define P_MUL(x, a, b) lpoly_mul((a), (b))
#define P_ONE(x) kp_one()
#define P_ZERO(x) ((LPoly){ 0, 0, NULL })
#define P_FREE(p) lpoly_free(p)
#define P_COEF(p, k) (&(p)->a[(k) - (p)->lo])
#define P_LO(p) ((p)->lo)
#define P_N(p) ((p)->n)
#define C_MUL(x, a, b) qfn_mul((a), (b))
#define C_INV(x, a) qfn_inv(a)
#define C_FREE(a) qfn_free(a)
#define C_SCALE_INT(x, a, k) qfn_scale_int((a), (k))
#define P_TAKE(x, arr, lo, n) lpoly_take((arr), (lo), (n))
#define C_ADD(x, a, b, sign) qfn_add((a), (b), (sign))
#define C_ZERO(x) qfn_zero()
#define C_IS_ZERO(a) qfn_is_zero(a)
#define C_COPY(a) qfn_copy(a)
#include "kpoly_tmpl.h"

LPoly kp_copy(const LPoly* a) { return kq_copy(0, a); }

LPoly kp_shift(const LPoly* a, int k) {
    LPoly r = kp_copy(a);
    if (r.n) r.lo += k;
    return r;
}

LPoly kp_scale(const LPoly* a, const QFn* c) { return kq_scale(0, a, c); }
LPoly kp_monic(const LPoly* a) { return kq_monic(0, a); }
LPoly kp_dtheta(const LPoly* a) { return kq_dtheta(0, a); }
KDiv kp_divmod(const LPoly* a, const LPoly* b) { return kq_divmod(0, a, b); }
LPoly kp_exact_div(const LPoly* a, const LPoly* b) { return kq_exact_div(0, a, b); }
LPoly kp_gcd(const LPoly* a, const LPoly* b) { return kq_gcd(0, a, b); }
KBezout kp_bezout(const LPoly* a, const LPoly* b, const LPoly* c) { return kq_bezout(0, a, b, c); }
KSquareFree kp_square_free(const LPoly* d) { return kq_square_free(0, d); }

void kdiv_free(KDiv* d) { lpoly_free(&d->q); lpoly_free(&d->r); }
void kbezout_free(KBezout* z) { lpoly_free(&z->s); lpoly_free(&z->t); }

void ksf_free(KSquareFree* s) {
    for (int i = 0; i < s->m; i++) lpoly_free(&s->f[i]);
    free(s->f);
    s->f = NULL;
    s->m = 0;
}

// --- fractions n/d in K(θ) ---

KFrac kfrac_make(LPoly n, LPoly d, TowerKind kind) {
    if (kind == TOWER_EXP && d.n) { n.lo -= d.lo; d.lo = 0; } // θ is a unit in K[θ, 1/θ]
    int s = kind == TOWER_EXP && n.n ? n.lo : 0;          // likewise for the numerator's θ^lo
    LPoly n0 = kp_shift(&n, -s), g = kp_gcd(&n0, &d);
    LPoly nq = kp_exact_div(&n0, &g), dq = kp_exact_div(&d, &g);
    QFn ic = qfn_inv(&dq.a[dq.n - 1]);                    // monic denominator
    LPoly ns = kp_scale(&nq, &ic), ds = kp_scale(&dq, &ic);
    KFrac f = { kp_shift(&ns, s), ds };
    lpoly_free(&n); lpoly_free(&d); lpoly_free(&n0); lpoly_free(&g);
    lpoly_free(&nq); lpoly_free(&dq); lpoly_free(&ns);
    qfn_free(&ic);
    return f;
}

void kfrac_free(KFrac* f) { lpoly_free(&f->n); lpoly_free(&f->d); }

// D(n/d) = (D(n)·d - n·D(d)) / d²
KFrac kfrac_derivative(const KFrac* f, const Tower* t) {
    LPoly dn = lpoly_derivative(&f->n, t), dd = lpoly_derivative(&f->d, t);
    LPoly l = lpoly_mul(&dn, &f->d), r = lpoly_mul(&f->n, &dd);
    KFrac out = kfrac_make(lpoly_add(&l, &r, -1), lpoly_mul(&f->d, &f->d), t->kind);
    lpoly_free(&dn); lpoly_free(&dd); lpoly_free(&l); lpoly_free(&r);
    return out;
}

// Henrici's formulas: with reduced inputs, only gcds of small pieces are needed (one big gcd
// of the result would bring Euclid's coefficient growth over Q(x) back)

// exp's numerators are Laurent: divisibility by a θ-free d ignores the θ^lo factor; for ln, θ is
// a normal polynomial like any other and stays
static int theta_shift(const LPoly* n, TowerKind kind) { return kind == TOWER_EXP && n->n ? n->lo : 0; }

// gcd(n, d) for a numerator n and a denominator d
static LPoly gcd_num(const LPoly* n, const LPoly* d, TowerKind kind) {
    LPoly n0 = kp_shift(n, -theta_shift(n, kind)), g = kp_gcd(&n0, d);
    lpoly_free(&n0);
    return g;
}

// n/g for g | n (up to exp's θ^lo)
static LPoly div_num(const LPoly* n, const LPoly* g, TowerKind kind) {
    int s = theta_shift(n, kind);
    LPoly n0 = kp_shift(n, -s), q = kp_exact_div(&n0, g), r = kp_shift(&q, s);
    lpoly_free(&n0); lpoly_free(&q);
    return r;
}

static KFrac kfrac_zero(void) { return (KFrac){ { 0, 0, NULL }, kp_one() }; }

typedef struct {
    const KFrac *a, *b;
    int sign;
} KSum;

// a/b ± c/d = (a·d' ± c·b')/(b·d') with g = gcd(b, d), b' = b/g, d' = d/g; the only common
// factor left can come from g
static KFrac kfrac_sum(const KSum* p, TowerKind kind) {
    const KFrac *a = p->a, *b = p->b;
    LPoly g = kp_gcd(&a->d, &b->d), a1 = kp_exact_div(&a->d, &g), b1 = kp_exact_div(&b->d, &g);
    LPoly l = lpoly_mul(&a->n, &b1), r = lpoly_mul(&b->n, &a1), num = lpoly_add(&l, &r, p->sign);
    KFrac out = kfrac_zero();
    if (num.n) {
        LPoly den = lpoly_mul(&a->d, &b1), h = gcd_num(&num, &g, kind);
        lpoly_free(&out.d);
        out = (KFrac){ div_num(&num, &h, kind), kp_exact_div(&den, &h) };
        lpoly_free(&den); lpoly_free(&h);
    }
    lpoly_free(&g); lpoly_free(&a1); lpoly_free(&b1); lpoly_free(&l); lpoly_free(&r); lpoly_free(&num);
    return out;
}

// (a/b)·(c/d) = (a/g1)·(c/g2) / ((b/g2)·(d/g1)) with g1 = gcd(a, d), g2 = gcd(c, b)
KFrac kfrac_mul(const KFrac* x, const KFrac* y, TowerKind kind) {
    if (!x->n.n || !y->n.n) return kfrac_zero();
    LPoly g1 = gcd_num(&x->n, &y->d, kind), g2 = gcd_num(&y->n, &x->d, kind);
    LPoly n1 = div_num(&x->n, &g1, kind), n2 = div_num(&y->n, &g2, kind);
    LPoly d1 = kp_exact_div(&x->d, &g2), d2 = kp_exact_div(&y->d, &g1);
    KFrac out = { lpoly_mul(&n1, &n2), lpoly_mul(&d1, &d2) };
    lpoly_free(&g1); lpoly_free(&g2); lpoly_free(&n1); lpoly_free(&n2); lpoly_free(&d1); lpoly_free(&d2);
    return out;
}

// d/n, already coprime: only the normalization (θ powers out of the denominator, monic)
KFrac kfrac_inv(const KFrac* f, TowerKind kind) {
    int s = theta_shift(&f->n, kind);
    LPoly n = kp_shift(&f->d, -s), d = kp_shift(&f->n, -s);
    QFn ic = qfn_inv(&d.a[d.n - 1]);
    KFrac out = { kp_scale(&n, &ic), kp_scale(&d, &ic) };
    lpoly_free(&n); lpoly_free(&d); qfn_free(&ic);
    return out;
}

KFrac kfrac_add(const KFrac* a, const KFrac* b, TowerKind kind) { return kfrac_sum(&(KSum){ a, b, 1 }, kind); }
KFrac kfrac_sub(const KFrac* a, const KFrac* b, TowerKind kind) { return kfrac_sum(&(KSum){ a, b, -1 }, kind); }

// --- integer forms and printing ---

LPoly kp_pow(const LPoly* a, int k) {
    LPoly r = kp_one();
    for (int i = 0; i < k; i++) { LPoly p = lpoly_mul(&r, a); lpoly_free(&r); r = p; }
    return r;
}

// gcd of two nonnegative rationals: gcd of the numerators over lcm of the denominators
static void rat_gcd(Rational* g, const Rational* a) {
    if (rat_is_zero(g)) { rat_set(g, a); return; }
    mpz_t n, d;
    mpz_inits(n, d, NULL);
    mpz_gcd(n, mpq_numref(g->val), mpq_numref(a->val));
    mpz_lcm(d, mpq_denref(g->val), mpq_denref(a->val));
    mpq_set_num(g->val, n);
    mpq_set_den(g->val, d);
    mpq_canonicalize(g->val);
    mpz_clears(n, d, NULL);
}

// lcm of the denominators of a's and b's coefficients (a monic polynomial in x); b may be NULL
static Polynomial* den_lcm(const LPoly* a, const LPoly* b) {
    Rational one = rat_one_val();
    Polynomial* l = poly_constant(&one);
    rat_clear(&one);
    const LPoly* both[] = { a, b };
    for (int j = 0; j < 2 && both[j]; j++)
        for (int i = 0; i < both[j]->n; i++) { Polynomial* m = poly_lcm(l, both[j]->a[i].den); free_polynomial(l); l = m; }
    return l;
}

// c·a for a polynomial c in x
static LPoly kp_times_x(const LPoly* a, Polynomial* c) {
    Rational one = rat_one_val();
    QFn q = { c, poly_constant(&one) };
    LPoly r = kp_scale(a, &q);
    qfn_free(&q);
    rat_clear(&one);
    return r;
}

LPoly kp_primitive(const LPoly* a) {
    LPoly b = kp_times_x(a, den_lcm(a, NULL));            // coefficients in Q[x]
    Rational g, c;
    rat_init(&g);
    for (int i = 0; i < b.n; i++) { c = poly_content(b.a[i].num); rat_gcd(&g, &c); rat_clear(&c); }
    Rational lc;
    rat_init(&lc);
    const Polynomial* top = b.a[b.n - 1].num;
    poly_coeff(&lc, top, polynomial_degree(top));
    if (rat_is_negative(&lc)) rat_neg(&g, &g);           // positive leading coefficient
    mpq_inv(g.val, g.val);
    LPoly r = kp_times_x(&b, poly_constant(&g));
    lpoly_free(&b);
    rat_clear(&g); rat_clear(&lc);
    return r;
}

// n/d with the x-denominators cleared: (Σ n_k·θ^k) / (Σ d_k·θ^k), n_k, d_k in Q[x]
Expr* kfrac_expr(const KFrac* f, const Tower* t) {
    if (kp_is_one(&f->d)) return lpoly_expr(&f->n, t);
    Polynomial* l = den_lcm(&f->n, &f->d);
    LPoly n = kp_times_x(&f->n, poly_copy(l)), d = kp_times_x(&f->d, l);
    Expr* e = create_div(lpoly_expr_order(&n, t, n.lo >= 0), lpoly_expr_order(&d, t, true));
    lpoly_free(&n); lpoly_free(&d);
    return e;
}

bool qfn_constant(const QFn* a, Rational* c) {
    if (polynomial_degree(a->num) > 0 || polynomial_degree(a->den) != 0) return false;
    poly_coeff(c, a->num, 0);                             // the denominator is monic: 1
    return true;
}

QFn qfn_rat(const Rational* c) {
    Rational one = rat_one_val();
    QFn q = { poly_constant(c), poly_constant(&one) };
    rat_clear(&one);
    return q;
}

QFn qfn_const(long n, long d) {
    Rational c;
    rat_init(&c);
    rat_set_si(&c, n, d);
    QFn q = qfn_rat(&c);
    rat_clear(&c);
    return q;
}

QFn kp_coeff(const LPoly* f, int k) {
    int i = k - f->lo;
    return i >= 0 && i < f->n ? qfn_copy(&f->a[i]) : qfn_zero();
}

KPair kp_common_multiple(const LPoly* a, const LPoly* b) {
    Polynomial* l = den_lcm(a, b);
    KPair r = { kp_times_x(a, poly_copy(l)), kp_times_x(b, l) };
    return r;
}
