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

LPoly kp_copy(const LPoly* a) {
    LPoly z = { 0, 0, NULL };
    return lpoly_add(a, &z, 1);
}

LPoly kp_shift(const LPoly* a, int k) {
    LPoly r = kp_copy(a);
    if (r.n) r.lo += k;
    return r;
}

LPoly kp_scale(const LPoly* a, const QFn* c) {
    LPoly t = lpoly_term(qfn_copy(c), 0), r = lpoly_mul(a, &t);
    lpoly_free(&t);
    return r;
}

static const QFn* kp_lc(const LPoly* a) { return &a->a[a->n - 1]; }

LPoly kp_monic(const LPoly* a) {
    QFn inv = qfn_inv(kp_lc(a));
    LPoly r = kp_scale(a, &inv);
    qfn_free(&inv);
    return r;
}

LPoly kp_dtheta(const LPoly* a) {
    LPoly r = { 0, 0, NULL };
    for (int i = 0; i < a->n; i++) {
        int k = a->lo + i;
        if (k == 0) continue;
        Rational kk = rat_from_int(k);
        QFn c = { poly_scale(a->a[i].num, &kk), poly_copy(a->a[i].den) };
        LPoly t = lpoly_term(c, k - 1), s = lpoly_add(&r, &t, 1);
        lpoly_free(&r); lpoly_free(&t);
        r = s;
        rat_clear(&kk);
    }
    return r;
}

// a = q·b + r with deg r < deg b, from the top one term at a time
KDiv kp_divmod(const LPoly* a, const LPoly* b) {
    KDiv d = { { 0, 0, NULL }, kp_copy(a) };
    QFn ib = qfn_inv(kp_lc(b));
    while (kp_deg(&d.r) >= kp_deg(b)) {
        LPoly t = lpoly_term(qfn_mul(kp_lc(&d.r), &ib), kp_deg(&d.r) - kp_deg(b)); // next quotient term
        LPoly tb = lpoly_mul(&t, b), r = lpoly_add(&d.r, &tb, -1), q = lpoly_add(&d.q, &t, 1);
        lpoly_free(&d.r); lpoly_free(&d.q); lpoly_free(&t); lpoly_free(&tb);
        d.r = r; d.q = q;
    }
    qfn_free(&ib);
    return d;
}

LPoly kp_exact_div(const LPoly* a, const LPoly* b) {
    KDiv d = kp_divmod(a, b);
    lpoly_free(&d.r);
    return d.q;
}

LPoly kp_gcd(const LPoly* a, const LPoly* b) {
    LPoly r0 = kp_copy(a), r1 = kp_copy(b);
    while (r1.n) {                                        // (r0, r1) <- (r1, r0 mod r1)
        KDiv d = kp_divmod(&r0, &r1);
        lpoly_free(&d.q); lpoly_free(&r0);
        r0 = r1; r1 = d.r;
    }
    LPoly g = r0.n ? kp_monic(&r0) : r0;
    if (r0.n) lpoly_free(&r0);
    lpoly_free(&r1);                                      // the zero remainder still owns its array
    return g;
}

// s with s·a ≡ 1 (mod b), by the extended Euclidean algorithm (only the s column is needed)
static LPoly kp_inverse(const LPoly* a, const LPoly* b) {
    LPoly r0 = kp_copy(a), r1 = kp_copy(b), s0 = kp_one(), s1 = { 0, 0, NULL };
    while (r1.n) {
        KDiv d = kp_divmod(&r0, &r1);
        LPoly qs = lpoly_mul(&d.q, &s1), s2 = lpoly_add(&s0, &qs, -1);
        lpoly_free(&qs); lpoly_free(&d.q); lpoly_free(&r0); lpoly_free(&s0);
        r0 = r1; r1 = d.r;
        s0 = s1; s1 = s2;
    }
    QFn inv = qfn_inv(&r0.a[0]);                          // r0 is the constant gcd
    LPoly s = kp_scale(&s0, &inv);
    qfn_free(&inv);
    lpoly_free(&r0); lpoly_free(&r1); lpoly_free(&s0); lpoly_free(&s1);
    return s;
}

KBezout kp_bezout(const LPoly* a, const LPoly* b, const LPoly* c) {
    LPoly inv = kp_inverse(a, b), cs = lpoly_mul(c, &inv);
    KDiv sd = kp_divmod(&cs, b);                          // s = c·a⁻¹ mod b
    LPoly sa = lpoly_mul(&sd.r, a), rest = lpoly_add(c, &sa, -1);
    KBezout z = { sd.r, kp_exact_div(&rest, b) };         // t = (c - s·a)/b
    lpoly_free(&inv); lpoly_free(&cs); lpoly_free(&sd.q); lpoly_free(&sa); lpoly_free(&rest);
    return z;
}

// Yun: with c = gcd(d, d_θ), w = d/c, y = d_θ/c, each step peels off the part of multiplicity i
KSquareFree kp_square_free(const LPoly* d) {
    KSquareFree s = { malloc((size_t)(kp_deg(d) + 1) * sizeof(LPoly)), 0 };
    LPoly dp = kp_dtheta(d), c = kp_gcd(d, &dp);
    LPoly w = kp_exact_div(d, &c), y = kp_exact_div(&dp, &c);
    lpoly_free(&dp); lpoly_free(&c);
    while (kp_deg(&w) > 0) {
        LPoly wp = kp_dtheta(&w), z = lpoly_add(&y, &wp, -1);
        LPoly g = kp_gcd(&w, &z);
        LPoly nw = kp_exact_div(&w, &g), ny = kp_exact_div(&z, &g);
        s.f[s.m++] = g;                                    // the part of multiplicity m
        lpoly_free(&wp); lpoly_free(&z); lpoly_free(&w); lpoly_free(&y);
        w = nw; y = ny;
    }
    lpoly_free(&w); lpoly_free(&y);
    return s;
}

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

typedef struct {
    const KFrac *a, *b;
    int sign;
} KSum;

// a ± b
static KFrac kfrac_sum(const KSum* p, TowerKind kind) {
    const KFrac *a = p->a, *b = p->b;
    int sign = p->sign;
    LPoly l = lpoly_mul(&a->n, &b->d), r = lpoly_mul(&b->n, &a->d);
    KFrac out = kfrac_make(lpoly_add(&l, &r, sign), lpoly_mul(&a->d, &b->d), kind);
    lpoly_free(&l); lpoly_free(&r);
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
