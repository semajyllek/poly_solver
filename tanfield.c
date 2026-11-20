// Conversion of rational functions of x, sin and cos into Q(x)(tan(H)). Elements are first
// written as a + cos(H)·b with a, b ∈ Q(x)(θ) (cos(H) itself is not rational in θ, only its
// square is); f lies in Q(x)(θ) exactly when b = 0.
#include "tanfield.h"
#include "expsum.h"
#include <stdlib.h>

#define MAX_MULTIPLE 64 // |j| in sin(j·H)

// a + cos(H)·b
typedef struct {
    KFrac a, b;
} TPair;

typedef struct {
    const Polynomial* H;
    const char* unsupported;
} TConv;

static KFrac kf_zero(void) { return (KFrac){ { 0, 0, NULL }, kp_one() }; }
static KFrac kf_copy(const KFrac* a) { return (KFrac){ kp_copy(&a->n), kp_copy(&a->d) }; }
static bool kf_is_zero(const KFrac* a) { return a->n.n == 0; }

static TPair tp_zero(void) { return (TPair){ kf_zero(), kf_zero() }; }

// -a (still reduced)
static KFrac kf_neg(const KFrac* a) {
    QFn m1 = qfn_const(-1, 1);
    KFrac r = { kp_scale(&a->n, &m1), kp_copy(&a->d) };
    qfn_free(&m1);
    return r;
}

static void tp_free(TPair* p) { kfrac_free(&p->a); kfrac_free(&p->b); }

LPoly tan_special(void) {
    LPoly one = kp_one(), t2 = kp_shift(&one, 2), r = lpoly_add(&one, &t2, 1);
    lpoly_free(&one); lpoly_free(&t2);
    return r;
}

// cos(H)² = 1/(1 + θ²)
static KFrac cos_squared(void) { return kfrac_make(kp_one(), tan_special(), TOWER_TAN); }

static TPair tp_add(const TPair* x, const TPair* y, int sign) {
    if (sign > 0) return (TPair){ kfrac_add(&x->a, &y->a, TOWER_TAN), kfrac_add(&x->b, &y->b, TOWER_TAN) };
    return (TPair){ kfrac_sub(&x->a, &y->a, TOWER_TAN), kfrac_sub(&x->b, &y->b, TOWER_TAN) };
}

// u·v, with zero factors skipped (b = 0 throughout for θ = tan(B/2))
static KFrac kf_mul(const KFrac* u, const KFrac* v) { return kf_is_zero(u) || kf_is_zero(v) ? kf_zero() : kfrac_mul(u, v, TOWER_TAN); }

// x.a·y.a + c²·x.b·y.b, the part of x·y without cos(H)
static KFrac even_product(const TPair* x, const TPair* y) {
    if (kf_is_zero(&x->b) || kf_is_zero(&y->b)) return kf_mul(&x->a, &y->a);
    KFrac c2 = cos_squared(), uv = kf_mul(&x->a, &y->a), wz = kfrac_mul(&x->b, &y->b, TOWER_TAN);
    KFrac cwz = kfrac_mul(&c2, &wz, TOWER_TAN), r = kfrac_add(&uv, &cwz, TOWER_TAN);
    kfrac_free(&c2); kfrac_free(&uv); kfrac_free(&wz); kfrac_free(&cwz);
    return r;
}

// (a + c·b)(a' + c·b') = (a·a' + c²·b·b') + c·(a·b' + b·a')
static TPair tp_mul(const TPair* x, const TPair* y) {
    KFrac ab = kf_mul(&x->a, &y->b), ba = kf_mul(&x->b, &y->a);
    TPair r = { even_product(x, y), kfrac_add(&ab, &ba, TOWER_TAN) };
    kfrac_free(&ab); kfrac_free(&ba);
    return r;
}

// 1/(a + c·b) = (a - c·b)/(a² - c²·b²)
static TPair tp_inv(const TPair* x, TConv* C) {
    KFrac nb = kf_neg(&x->b);
    KFrac n = even_product(x, &(TPair){ x->a, nb });      // (a + c·b)(a - c·b), read only
    TPair r;
    if (kf_is_zero(&n)) { C->unsupported = "division by zero"; r = tp_zero(); }
    else {
        KFrac in = kfrac_inv(&n, TOWER_TAN);
        r = (TPair){ kf_mul(&x->a, &in), kf_mul(&nb, &in) };
        kfrac_free(&in);
    }
    kfrac_free(&nb); kfrac_free(&n);
    return r;
}

// (1 + iθ)^n = re + i·im
static void binomial_parts(int n, LPoly* re, LPoly* im) {
    *re = (LPoly){ 0, 0, NULL };
    *im = (LPoly){ 0, 0, NULL };
    for (int k = 0; k <= n; k++) {
        Rational c;
        rat_init(&c);
        mpz_bin_uiui(mpq_numref(c.val), (unsigned long)n, (unsigned long)k);
        if (k % 4 >= 2) rat_neg(&c, &c);                  // i^k = 1, i, -1, -i
        LPoly t = lpoly_term(qfn_rat(&c), k), *acc = k % 2 ? im : re, s = lpoly_add(acc, &t, 1);
        lpoly_free(acc); lpoly_free(&t);
        *acc = s;
        rat_clear(&c);
    }
}

// sin(jH) or cos(jH) = Im or Re of (1 + iθ)^|j|, times cos(H)^|j| (and -1 for sin with j < 0)
static TPair leaf(long j, bool sine) {
    int n = (int)labs(j);
    LPoly re, im, den = kp_one(), q = tan_special();
    binomial_parts(n, &re, &im);
    for (int i = 0; i < n / 2; i++) { LPoly t = lpoly_mul(&den, &q); lpoly_free(&den); den = t; }
    QFn sign = qfn_const(sine && j < 0 ? -1 : 1, 1);
    KFrac v = kfrac_make(kp_scale(sine ? &im : &re, &sign), den, TOWER_TAN);
    lpoly_free(&re); lpoly_free(&im); lpoly_free(&q); qfn_free(&sign);
    return n % 2 ? (TPair){ kf_zero(), v } : (TPair){ v, kf_zero() };
}

// the polynomial a sin/cos argument must be
static Polynomial* trig_arg(const Expr* arg) {
    Expr* c = canonicalize((Expr*)arg);                   // canonicalize does not modify arg
    Polynomial* p = c && c->type == EXPR_POLY && polynomial_degree(c->poly) > 0 ? poly_copy(c->poly) : NULL;
    free_expr(c);
    return p;
}

// lc(p)/lc(q)
static void lc_ratio(Rational* r, const Polynomial* p, const Polynomial* q) {
    Rational b;
    rat_init(&b);
    poly_coeff(r, p, polynomial_degree(p));
    poly_coeff(&b, q, polynomial_degree(q));
    rat_div(r, r, &b);
    rat_clear(&b);
}

// p = r·q for a rational r
static bool proportional(const Polynomial* p, const Polynomial* q, Rational* r) {
    lc_ratio(r, p, q);
    Polynomial *rq = poly_scale(q, r), *d = poly_subtract(rq, p);
    bool same = poly_is_zero(d);
    free_polynomial(rq); free_polynomial(d);
    return same;
}

static TPair conv(const Expr* e, TConv* C);

static TPair conv_trig(const Expr* e, TConv* C) {
    Polynomial* u = trig_arg(e->fn_arg);
    Rational r;
    rat_init(&r);
    bool ok = u && proportional(u, C->H, &r) && mpz_cmp_ui(mpq_denref(r.val), 1) == 0
           && mpz_cmpabs_ui(mpq_numref(r.val), MAX_MULTIPLE) <= 0;
    long j = ok ? mpz_get_si(mpq_numref(r.val)) : 0;
    free_polynomial(u);
    rat_clear(&r);
    if (!ok) { C->unsupported = "sin and cos here must take integer multiples of one polynomial (such as sin(x) and cos(2x))"; return tp_zero(); }
    return leaf(j, e->type == EXPR_SIN);
}

static TPair conv_pow(const Expr* e, TConv* C) {
    long k = 0;
    if (!small_integer(e->binop.right, &k)) { C->unsupported = "exponents must be small integers"; return tp_zero(); }
    TPair b = conv(e->binop.left, C);
    if (k < 0 && !C->unsupported) { TPair inv = tp_inv(&b, C); tp_free(&b); b = inv; k = -k; }
    TPair r = { kfrac_make(kp_one(), kp_one(), TOWER_TAN), kf_zero() };
    for (long i = 0; i < k && !C->unsupported; i++) { TPair p = tp_mul(&r, &b); tp_free(&r); r = p; }
    tp_free(&b);
    return r;
}

static TPair conv_binop(const Expr* e, TConv* C) {
    if (e->binop.op == OP_POW) return conv_pow(e, C);
    TPair l = conv(e->binop.left, C), r = conv(e->binop.right, C), out = tp_zero();
    if (!C->unsupported) {
        tp_free(&out);
        switch (e->binop.op) {
            case OP_ADD: out = tp_add(&l, &r, 1); break;
            case OP_SUB: out = tp_add(&l, &r, -1); break;
            case OP_MUL: out = tp_mul(&l, &r); break;
            default: { TPair inv = tp_inv(&r, C); out = tp_mul(&l, &inv); tp_free(&inv); }
        }
    }
    tp_free(&l); tp_free(&r);
    return out;
}

// a rational function of x
static TPair conv_rational(const Expr* e, TConv* C) {
    Expr* c = canonicalize((Expr*)e);
    if (!c) { C->unsupported = "not a rational function of x, sin and cos"; return tp_zero(); }
    QFn a = c->type == EXPR_POLY ? qfn_poly(c->poly) : qfn_from(poly_copy(c->ratfn.num), poly_copy(c->ratfn.den));
    free_expr(c);
    return (TPair){ kfrac_make(lpoly_term(a, 0), kp_one(), TOWER_TAN), kf_zero() };
}

static TPair conv(const Expr* e, TConv* C) {
    if (C->unsupported) return tp_zero();
    if (expr_has_exp(e) || expr_has_ln(e)) {
        C->unsupported = "exp or ln together with a sin/cos denominator is not supported yet";
        return tp_zero();
    }
    if (!expr_has_trig(e)) return conv_rational(e, C);
    switch (e->type) {
        case EXPR_BINOP: return conv_binop(e, C);
        case EXPR_NEG: { TPair f = conv(e->operand, C), r = { kf_neg(&f.a), kf_neg(&f.b) }; tp_free(&f); return r; }
        case EXPR_SIN: case EXPR_COS: return conv_trig(e, C);
        default: C->unsupported = "only +, -, *, /, integer powers, sin and cos are supported here"; return tp_zero();
    }
}

// the first sin/cos argument (a nonconstant polynomial), or NULL
static Polynomial* first_arg(const Expr* e) {
    switch (e->type) {
        case EXPR_SIN: case EXPR_COS: return trig_arg(e->fn_arg);
        case EXPR_NEG: return first_arg(e->operand);
        case EXPR_BINOP: {
            Polynomial* p = first_arg(e->binop.left);
            return p ? p : first_arg(e->binop.right);
        }
        default: return NULL;
    }
}

// the rational gcd g of the multiples r (A = r·A0 for each argument A), folded in
static void fold_multiple(const Expr* e, const Polynomial* A0, Rational* g, bool* ok) {
    if (!*ok) return;
    if (e->type == EXPR_NEG) { fold_multiple(e->operand, A0, g, ok); return; }
    if (e->type == EXPR_BINOP) { fold_multiple(e->binop.left, A0, g, ok); fold_multiple(e->binop.right, A0, g, ok); return; }
    if (e->type != EXPR_SIN && e->type != EXPR_COS) return;
    Polynomial* u = trig_arg(e->fn_arg);
    Rational r;
    rat_init(&r);
    *ok = u && proportional(u, A0, &r);
    if (*ok) {                                             // gcd(a/b, c/d) = gcd(a, c)/lcm(b, d)
        mpq_abs(r.val, r.val);
        mpz_gcd(mpq_numref(g->val), mpq_numref(g->val), mpq_numref(r.val));
        mpz_lcm(mpq_denref(g->val), mpq_denref(g->val), mpq_denref(r.val));
        mpq_canonicalize(g->val);
    }
    free_polynomial(u);
    rat_clear(&r);
}

// B: every sin/cos argument is an integer multiple of it; NULL when there is none
static Polynomial* base_of(const Expr* e) {
    Polynomial* A0 = first_arg(e);
    if (!A0) return NULL;
    Rational g;
    rat_init(&g);                                          // 0/1: gcd(0, r) = r
    bool ok = true;
    fold_multiple(e, A0, &g, &ok);
    Rational lc;
    rat_init(&lc);
    poly_coeff(&lc, A0, polynomial_degree(A0));
    if (rat_is_negative(&lc)) rat_neg(&g, &g);             // lc(B) > 0: tan(x/2), not tan(-x/2)
    rat_clear(&lc);
    Polynomial* B = ok ? poly_scale(A0, &g) : NULL;
    free_polynomial(A0);
    rat_clear(&g);
    return B;
}

// the parity of e under B -> B + π (sin(kB), cos(kB) have k mod 2), or -1 when a sum mixes
// parities; parity 0 means e is rational in tan(B) (H = B base, multiples k known by C)
static int parity(const Expr* e, const TConv* C) {
    if (!expr_has_trig(e)) return 0;
    switch (e->type) {
        case EXPR_NEG: return parity(e->operand, C);
        case EXPR_SIN: case EXPR_COS: {
            Polynomial* u = trig_arg(e->fn_arg);
            Rational r;
            rat_init(&r);
            int p = u && proportional(u, C->H, &r) && mpz_cmp_ui(mpq_denref(r.val), 1) == 0 ? (int)mpz_odd_p(mpq_numref(r.val)) : -1;
            free_polynomial(u);
            rat_clear(&r);
            return p;
        }
        case EXPR_BINOP: {
            int l = parity(e->binop.left, C);
            if (e->binop.op == OP_POW) {
                long k = 0;
                return l < 0 || !small_integer(e->binop.right, &k) ? -1 : (int)(l & (k & 1));
            }
            int r = parity(e->binop.right, C);
            if (l < 0 || r < 0) return -1;
            if (e->binop.op == OP_ADD || e->binop.op == OP_SUB) return l == r ? l : -1;
            return l ^ r;                                 // products and quotients add parities
        }
        default: return -1;
    }
}

// e in θ = tan(H), or false (C->unsupported, or b ≠ 0: e is not rational in tan(H))
static bool convert(const Expr* e, TConv* C, KFrac* out) {
    TPair p = conv(e, C);
    bool ok = !C->unsupported && kf_is_zero(&p.b);
    if (ok) *out = kf_copy(&p.a);
    tp_free(&p);
    return ok;
}

ExpForm tan_form(const Expr* e) {
    ExpForm F = { { TOWER_TAN, NULL, NULL, NULL, NULL, NULL, false }, { 0, 0, NULL }, kp_one(), NULL };
    Polynomial* B = base_of(e);
    if (!B) { F.unsupported = "sin and cos here must take integer multiples of one polynomial (such as sin(x) and cos(2x))"; return F; }
    Rational half = rat_from_int(1);
    mpq_set_si(half.val, 1, 2);
    Polynomial* H2 = poly_scale(B, &half);
    rat_clear(&half);
    TConv C = { B, NULL };
    KFrac f;
    if (parity(e, &C) != 0) C.H = H2;                     // odd or mixed: θ = tan(B/2), where b = 0 always
    bool ok = convert(e, &C, &f);                          // even (Bioche): θ = tan(B), and then b = 0 too
    if (ok) { lpoly_free(&F.den); F.f = f.n; F.den = f.d; F.t.v = poly_copy(C.H); F.t.half = C.H == H2; }
    else F.unsupported = C.unsupported ? C.unsupported : "internal error: not rational in tan(B/2)";
    free_polynomial(B); free_polynomial(H2);
    return F;
}

bool tan_frac(const Expr* e, const Polynomial* H, KFrac* f) {
    TConv C = { H, NULL };
    return convert(e, &C, f);
}

bool tan_poly(const Expr* e, const Polynomial* H, LPoly* p) {
    KFrac f;
    if (!tan_frac(e, H, &f)) return false;
    bool ok = kp_is_one(&f.d);
    if (ok) *p = kp_copy(&f.n);
    kfrac_free(&f);
    return ok;
}
