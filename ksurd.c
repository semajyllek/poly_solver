// K(√M)[θ]: just enough for a monic gcd (remainders and the leading-coefficient inverse).
#include "ksurd.h"

SPoly sp_from(const LPoly* p, const LPoly* q) { return (SPoly){ kp_copy(p), kp_copy(q) }; }

void sp_free(SPoly* a) { lpoly_free(&a->p); lpoly_free(&a->q); }

int sp_deg(const SPoly* a) {
    int dp = kp_deg(&a->p), dq = kp_deg(&a->q);
    return dp > dq ? dp : dq;
}

static void se_free(SElt* x) { qfn_free(&x->u); qfn_free(&x->w); }

// (u1 + √M·w1)(u2 + √M·w2) = (u1·u2 + M·w1·w2) + √M·(u1·w2 + w1·u2)
static SElt se_mul(const SElt* x, const SElt* y, const Rational* M) {
    QFn uu = qfn_mul(&x->u, &y->u), ww = qfn_mul(&x->w, &y->w), uw = qfn_mul(&x->u, &y->w), wu = qfn_mul(&x->w, &y->u);
    QFn m = qfn_rat(M), mww = qfn_mul(&m, &ww);
    SElt r = { qfn_add(&uu, &mww, 1), qfn_add(&uw, &wu, 1) };
    qfn_free(&uu); qfn_free(&ww); qfn_free(&uw); qfn_free(&wu); qfn_free(&m); qfn_free(&mww);
    return r;
}

// 1/(u + √M·w) = (u - √M·w)/(u² - M·w²); the norm is nonzero since M is not a square
static SElt se_inv(const SElt* x, const Rational* M) {
    QFn uu = qfn_mul(&x->u, &x->u), ww = qfn_mul(&x->w, &x->w), m = qfn_rat(M), mww = qfn_mul(&m, &ww);
    QFn n = qfn_add(&uu, &mww, -1), in = qfn_inv(&n), z = qfn_zero(), nw = qfn_add(&z, &x->w, -1);
    SElt r = { qfn_mul(&x->u, &in), qfn_mul(&nw, &in) };
    qfn_free(&uu); qfn_free(&ww); qfn_free(&m); qfn_free(&mww); qfn_free(&n); qfn_free(&in); qfn_free(&z); qfn_free(&nw);
    return r;
}

static SElt sp_coeff(const SPoly* a, int k) { return (SElt){ kp_coeff(&a->p, k), kp_coeff(&a->q, k) }; }

// c·θ^k
typedef struct {
    SElt c;
    int k;
} STerm;

// a·c·θ^k
static SPoly sp_mul_term(const SPoly* a, const STerm* t, const Rational* M) {
    LPoly cu = lpoly_term(qfn_copy(&t->c.u), t->k), cw = lpoly_term(qfn_copy(&t->c.w), t->k);
    LPoly pu = lpoly_mul(&a->p, &cu), qw = lpoly_mul(&a->q, &cw), pw = lpoly_mul(&a->p, &cw), qu = lpoly_mul(&a->q, &cu);
    QFn m = qfn_rat(M);
    LPoly mqw = kp_scale(&qw, &m);
    SPoly r = { lpoly_add(&pu, &mqw, 1), lpoly_add(&pw, &qu, 1) };
    lpoly_free(&cu); lpoly_free(&cw); lpoly_free(&pu); lpoly_free(&qw); lpoly_free(&pw); lpoly_free(&qu); lpoly_free(&mqw);
    qfn_free(&m);
    return r;
}

// a mod b, from the top one term at a time
static SPoly sp_rem(const SPoly* a, const SPoly* b, const Rational* M) {
    SPoly r = sp_from(&a->p, &a->q);
    int db = sp_deg(b);
    SElt lb = sp_coeff(b, db), ib = se_inv(&lb, M);
    while (sp_deg(&r) >= db) {
        int dr = sp_deg(&r);
        SElt lr = sp_coeff(&r, dr);
        STerm t = { se_mul(&lr, &ib, M), dr - db };
        SPoly tb = sp_mul_term(b, &t, M);
        SPoly n = { lpoly_add(&r.p, &tb.p, -1), lpoly_add(&r.q, &tb.q, -1) };
        sp_free(&r); sp_free(&tb); se_free(&lr); se_free(&t.c);
        r = n;
    }
    se_free(&lb); se_free(&ib);
    return r;
}

SPoly sp_gcd(const SPoly* a, const SPoly* b, const Rational* M) {
    SPoly r0 = sp_from(&a->p, &a->q), r1 = sp_from(&b->p, &b->q);
    while (sp_deg(&r1) >= 0) {                            // (r0, r1) <- (r1, r0 mod r1)
        SPoly r = sp_rem(&r0, &r1, M);
        sp_free(&r0);
        r0 = r1; r1 = r;
    }
    int d = sp_deg(&r0);
    SElt lc = sp_coeff(&r0, d), inv = se_inv(&lc, M);
    STerm t = { inv, 0 };
    SPoly g = sp_mul_term(&r0, &t, M);
    se_free(&lc); se_free(&inv);
    sp_free(&r0); sp_free(&r1);
    return g;
}
