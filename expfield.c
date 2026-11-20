// Q(x)(θ), θ = exp(v): rational functions, Laurent polynomials in θ, and the conversion
// from an expression tree (which also finds v).
#include "expfield.h"
#include "modgcd.h"
#include "kpoly.h"
#include <stdlib.h>

#define MAX_EXP_POWER 256 // |n| in u^n, like canonicalize's cap on polynomial powers

// --- Q(x) ---

static Polynomial* poly_one(void) {
    Rational one = rat_one_val();
    Polynomial* p = poly_constant(&one);
    rat_clear(&one);
    return p;
}

QFn qfn_zero(void) { return (QFn){ create_polynomial(), poly_one() }; }

QFn qfn_from(Polynomial* num, Polynomial* den) {
    if (poly_is_zero(num)) { free_polynomial(num); free_polynomial(den); return qfn_zero(); }
    Polynomial* g = poly_gcd_fast(num, den);
    Polynomial *n = poly_divide(num, g), *d = poly_divide(den, g);
    free_polynomial(num); free_polynomial(den); free_polynomial(g);
    Rational lc;                                          // monic denominator
    rat_init(&lc);
    poly_coeff(&lc, d, polynomial_degree(d));
    mpq_inv(lc.val, lc.val);
    QFn r = { poly_scale(n, &lc), poly_scale(d, &lc) };
    free_polynomial(n); free_polynomial(d);
    rat_clear(&lc);
    return r;
}

QFn qfn_copy(const QFn* a) { return (QFn){ poly_copy(a->num), poly_copy(a->den) }; }
void qfn_free(QFn* a) { free_polynomial(a->num); free_polynomial(a->den); a->num = a->den = NULL; }
bool qfn_is_zero(const QFn* a) { return poly_is_zero(a->num); }

QFn qfn_poly(const Polynomial* p) { return (QFn){ poly_copy(p), poly_one() }; }

// a ± b
QFn qfn_add(const QFn* a, const QFn* b, int sign) {
    Polynomial *l = poly_multiply(a->num, b->den), *r = poly_multiply(b->num, a->den);
    Polynomial* n = sign > 0 ? poly_add(l, r) : poly_subtract(l, r);
    free_polynomial(l); free_polynomial(r);
    return qfn_from(n, poly_multiply(a->den, b->den));
}

QFn qfn_mul(const QFn* a, const QFn* b) {
    return qfn_from(poly_multiply(a->num, b->num), poly_multiply(a->den, b->den));
}

// 1/a for a ≠ 0
QFn qfn_inv(const QFn* a) { return qfn_from(poly_copy(a->den), poly_copy(a->num)); }

// (n/d)' = (n'd - nd')/d²
QFn qfn_derivative(const QFn* a) {
    Polynomial *dn = poly_derivative(a->num), *dd = poly_derivative(a->den);
    Polynomial *l = poly_multiply(dn, a->den), *r = poly_multiply(a->num, dd);
    Polynomial* n = poly_subtract(l, r);
    free_polynomial(dn); free_polynomial(dd); free_polynomial(l); free_polynomial(r);
    return qfn_from(n, poly_multiply(a->den, a->den));
}

Expr* qfn_expr(const QFn* a) {
    if (polynomial_degree(a->den) == 0) return create_poly_expr(poly_copy(a->num));
    return create_div(create_poly_expr(poly_copy(a->num)), create_poly_expr(poly_copy(a->den)));
}



// --- Laurent polynomials in θ ---

static LPoly lp_empty(void) { return (LPoly){ 0, 0, NULL }; }

void lpoly_free(LPoly* f) {
    for (int i = 0; i < f->n; i++) qfn_free(&f->a[i]);
    free(f->a);
    *f = lp_empty();
}

// zero coefficients for powers lo .. lo + n - 1
static LPoly lp_span(int lo, int n) {
    LPoly f = { lo, n, malloc((size_t)(n > 0 ? n : 1) * sizeof(QFn)) };
    for (int i = 0; i < n; i++) f.a[i] = qfn_zero();
    return f;
}

// drop zero coefficients at both ends
static void lp_trim(LPoly* f) {
    int s = 0, e = f->n;
    while (s < e && qfn_is_zero(&f->a[s])) qfn_free(&f->a[s++]);
    while (e > s && qfn_is_zero(&f->a[e - 1])) qfn_free(&f->a[--e]);
    for (int i = s; i < e; i++) f->a[i - s] = f->a[i];
    f->lo += s;
    f->n = e - s;
    if (f->n == 0) f->lo = 0;
}

// a·θ^k; takes a
LPoly lpoly_term(QFn a, int k) {
    LPoly f = { k, 1, malloc(sizeof(QFn)) };
    f.a[0] = a;
    lp_trim(&f);
    return f;
}

LPoly lpoly_add(const LPoly* f, const LPoly* g, int sign) {
    if (f->n == 0 && g->n == 0) return lp_empty();
    int lo = f->n == 0 ? g->lo : g->n == 0 ? f->lo : (f->lo < g->lo ? f->lo : g->lo);
    int hi = f->n == 0 ? g->lo + g->n : g->n == 0 ? f->lo + f->n
           : (f->lo + f->n > g->lo + g->n ? f->lo + f->n : g->lo + g->n);
    LPoly r = lp_span(lo, hi - lo);
    for (int i = 0; i < f->n; i++) {
        QFn* t = &r.a[f->lo + i - lo];
        QFn s = qfn_add(t, &f->a[i], 1);
        qfn_free(t); *t = s;
    }
    for (int i = 0; i < g->n; i++) {
        QFn* t = &r.a[g->lo + i - lo];
        QFn s = qfn_add(t, &g->a[i], sign);
        qfn_free(t); *t = s;
    }
    lp_trim(&r);
    return r;
}

LPoly lpoly_mul(const LPoly* f, const LPoly* g) {
    if (f->n == 0 || g->n == 0) return lp_empty();
    LPoly r = lp_span(f->lo + g->lo, f->n + g->n - 1);
    for (int i = 0; i < f->n; i++)
        for (int j = 0; j < g->n; j++) {
            QFn p = qfn_mul(&f->a[i], &g->a[j]);
            QFn* t = &r.a[i + j];
            QFn s = qfn_add(t, &p, 1);
            qfn_free(t); qfn_free(&p); *t = s;
        }
    lp_trim(&r);
    return r;
}

// (a·θ^k)' = a'·θ^k + k·a·L·θ^(k-1) for θ = ln(u), L = u'/u
static LPoly log_derivative(const LPoly* f, const Tower* t) {
    LPoly r = lp_span(f->lo > 0 ? f->lo - 1 : 0, f->n + 1);
    QFn L = tower_L(t);
    for (int i = 0; i < f->n; i++) {
        int k = f->lo + i;
        QFn d = qfn_derivative(&f->a[i]);
        QFn* at = &r.a[k - r.lo];
        QFn s = qfn_add(at, &d, 1);
        qfn_free(at); *at = s; qfn_free(&d);
        if (k == 0) continue;
        Rational kk = rat_from_int(k);
        QFn kL = { poly_scale(L.num, &kk), poly_copy(L.den) };
        QFn p = qfn_mul(&kL, &f->a[i]);
        QFn* below = &r.a[k - 1 - r.lo];
        QFn s2 = qfn_add(below, &p, 1);
        qfn_free(below); *below = s2;
        qfn_free(&kL); qfn_free(&p);
        rat_clear(&kk);
    }
    qfn_free(&L);
    lp_trim(&r);
    return r;
}

// r += a·θ^k; takes a
static void lp_accumulate(LPoly* r, QFn a, int k) {
    LPoly term = lpoly_term(a, k), s = lpoly_add(r, &term, 1);
    lpoly_free(r); lpoly_free(&term);
    *r = s;
}

// (a·θ^k)' = a'·θ^k + k·v'·a·(θ^(k-1) + θ^(k+1)) for θ = tan(v)
static LPoly tan_derivative(const LPoly* f, const Tower* t) {
    LPoly r = lp_empty();
    QFn vp = { poly_derivative(t->v), poly_one() };
    for (int i = 0; i < f->n; i++) {
        int k = f->lo + i;
        lp_accumulate(&r, qfn_derivative(&f->a[i]), k);
        if (k == 0) continue;
        QFn kq = qfn_const(k, 1), kv = qfn_mul(&kq, &vp), c = qfn_mul(&kv, &f->a[i]);
        lp_accumulate(&r, qfn_copy(&c), k - 1);
        lp_accumulate(&r, c, k + 1);
        qfn_free(&kq); qfn_free(&kv);
    }
    qfn_free(&vp);
    return r;
}

LPoly lpoly_derivative(const LPoly* f, const Tower* t) {
    if (t->kind == TOWER_LOG) return log_derivative(f, t);
    if (t->kind == TOWER_TAN) return tan_derivative(f, t);
    LPoly r = lp_span(f->lo, f->n);
    QFn vp = { poly_derivative(t->v), poly_one() };
    for (int i = 0; i < f->n; i++) {
        Rational k = rat_from_int(f->lo + i);
        QFn d = qfn_derivative(&f->a[i]);
        Polynomial* kv = poly_scale(vp.num, &k);
        QFn kvq = { kv, poly_one() };
        QFn prod = qfn_mul(&kvq, &f->a[i]);
        qfn_free(&r.a[i]);
        r.a[i] = qfn_add(&d, &prod, 1);
        qfn_free(&d); qfn_free(&kvq); qfn_free(&prod);
        rat_clear(&k);
    }
    qfn_free(&vp);
    lp_trim(&r);
    return r;
}

bool lpoly_equal(const LPoly* f, const LPoly* g) {
    LPoly d = lpoly_add(f, g, -1);
    bool z = d.n == 0;
    lpoly_free(&d);
    return z;
}

QFn tower_L(const Tower* t) { return (QFn){ poly_copy(t->Ln), poly_copy(t->Ld) }; }

// exp(k·v), ln(u)^k, or tan(v)^k (as sin(v)/cos(v), which prints as tan(v))
static Expr* theta_power(const Tower* t, int k) {
    if (t->kind == TOWER_TAN) {
        Expr* tn = create_div(create_sin(create_poly_expr(poly_copy(t->v))), create_cos(create_poly_expr(poly_copy(t->v))));
        return k == 1 ? tn : create_pow(tn, create_const(k));
    }
    if (t->kind == TOWER_LOG) {
        QFn u = { t->un, t->ud };
        Expr* ln = create_ln(qfn_expr(&u));
        return k == 1 ? ln : create_pow(ln, create_const(k));
    }
    Rational kk = rat_from_int(k);
    Expr* e = create_exp(create_poly_expr(poly_scale(t->v, &kk)));
    rat_clear(&kk);
    return e;
}

// a constant polynomial c: -1 and 1 print without a factor
static bool constant_value(const Polynomial* p, Rational* c) {
    if (polynomial_degree(p) > 0) return false;
    poly_coeff(c, p, 0);
    return true;
}

// a·exp(k·v) printed as num·exp(k·v), (num)·exp(k·v), or (num·exp(k·v)) / den
Expr* lpoly_term_expr(const QFn* a, const Tower* t, int k) {
    if (k == 0) return qfn_expr(a);
    Expr* e = theta_power(t, k);
    Rational c;
    rat_init(&c);
    bool constant = constant_value(a->num, &c);
    bool minus_one = constant && mpz_cmp_si(mpq_numref(c.val), -1) == 0 && mpz_cmp_ui(mpq_denref(c.val), 1) == 0;
    if (!minus_one && !(constant && rat_is_one(&c))) e = create_mul(create_poly_expr(poly_copy(a->num)), e);
    rat_clear(&c);
    if (polynomial_degree(a->den) > 0) e = create_div(e, create_poly_expr(poly_copy(a->den)));
    return minus_one ? create_neg(e) : e;                 // outside, so sums print a - b
}

Expr* lpoly_expr(const LPoly* f, const Tower* t) { return lpoly_expr_order(f, t, t->kind != TOWER_EXP); }

// the terms from the top power down (ln(x)^3 first, exp(x) + 1), or from the bottom up
Expr* lpoly_expr_order(const LPoly* f, const Tower* t, bool down) {
    Expr* acc = NULL;
    for (int j = 0; j < f->n; j++) {
        int i = down ? f->n - 1 - j : j;
        if (qfn_is_zero(&f->a[i])) continue;
        Expr* term = lpoly_term_expr(&f->a[i], t, f->lo + i);
        acc = acc ? create_add(acc, term) : term;
    }
    return acc ? acc : create_const(0);
}

// --- the tower: θ = exp(v) with every exp argument an integer multiple of v ---

typedef struct {
    Polynomial** u; // the nonconstant exp arguments found
    int n;
    QFn* w;         // the nonconstant ln arguments found
    int nw;
    const char* unsupported;
} Args;

static void note(Args* A, const char* why) { if (!A->unsupported) A->unsupported = why; }

// an ln argument: a nonconstant rational function (ln(1) = 0 is fine; ln(c) is a new constant)
static void collect_ln(const Expr* e, Args* A) {
    Expr* c = canonicalize(e->ln_arg);
    if (!c) { note(A, "ln of something other than a rational function (such as ln(ln(x))) is not supported yet"); return; }
    QFn u = c->type == EXPR_POLY ? qfn_poly(c->poly) : qfn_from(poly_copy(c->ratfn.num), poly_copy(c->ratfn.den));
    free_expr(c);
    bool constant = polynomial_degree(u.num) <= 0 && polynomial_degree(u.den) == 0;
    Rational v;
    rat_init(&v);
    if (constant) poly_coeff(&v, u.num, 0);
    if (constant && !rat_is_one(&v)) note(A, "ln of a constant other than 1 (a transcendental number) is not supported yet");
    rat_clear(&v);
    if (constant) { qfn_free(&u); return; }
    A->w = realloc(A->w, (size_t)(A->nw + 1) * sizeof(QFn));
    A->w[A->nw++] = u;
}

static void collect(const Expr* e, Args* A) {
    switch (e->type) {
        case EXPR_BINOP: collect(e->binop.left, A); collect(e->binop.right, A); return;
        case EXPR_NEG: collect(e->operand, A); return;
        case EXPR_LN: collect_ln(e, A); return;
        case EXPR_EXP: break;
        default: return;
    }
    Expr* c = canonicalize(e->fn_arg);
    if (!c || c->type != EXPR_POLY) note(A, "exp of a non-polynomial argument (such as exp(1/x) or exp(exp(x))) is not supported yet");
    else if (polynomial_degree(c->poly) > 0) {
        A->u = realloc(A->u, (size_t)(A->n + 1) * sizeof(Polynomial*));
        A->u[A->n++] = poly_copy(c->poly);
    } else if (!poly_is_zero(c->poly)) note(A, "exp of a nonzero constant (a transcendental number) is not supported yet");
    free_expr(c);
}

static void lc_of(Rational* r, const Polynomial* p) { poly_coeff(r, p, polynomial_degree(p)); }

// v = w/l with l the lcm of the denominators of the ratios u_i/w (the ratio for w itself is 1,
// so every u_i becomes an integer multiple of v), oriented with lc(v) > 0; false when some
// u_i is not a rational multiple of w (several independent exponentials)
static bool make_tower(const Args* A, Tower* t) {
    const Polynomial* w = A->u[0];
    Rational lw, li, r;
    rat_init(&lw); rat_init(&li); rat_init(&r);
    lc_of(&lw, w);
    mpz_t l;
    mpz_init_set_ui(l, 1);
    bool ok = true;
    for (int i = 0; ok && i < A->n; i++) {
        lc_of(&li, A->u[i]);
        rat_div(&r, &li, &lw);                            // u_i = r·w exactly?
        Polynomial* rw = poly_scale(w, &r);
        Polynomial* d = poly_subtract(A->u[i], rw);
        ok = poly_is_zero(d);
        free_polynomial(rw); free_polynomial(d);
        mpz_lcm(l, l, mpq_denref(r.val));
    }
    if (ok) {
        mpq_set_ui(r.val, 1, 1);
        mpq_set_den(r.val, l);
        if (rat_is_negative(&lw)) rat_neg(&r, &r);        // lc(v) > 0
        t->v = poly_scale(w, &r);
    }
    mpz_clear(l);
    rat_clear(&lw); rat_clear(&li); rat_clear(&r);
    return ok;
}

// the k with u = k·v (u nonconstant and a multiple of v by construction)
static int power_of(const Tower* t, const Polynomial* u) {
    Rational a, b;
    rat_init(&a); rat_init(&b);
    lc_of(&a, u);
    lc_of(&b, t->v);
    rat_div(&a, &a, &b);
    int k = (int)mpz_get_si(mpq_numref(a.val));
    rat_clear(&a); rat_clear(&b);
    return k;
}

// --- the log tower: θ = ln(u), every ln argument an integer power of u ---

static QFn log_deriv(const QFn* u) {
    QFn d = qfn_derivative(u), inv = qfn_inv(u);
    QFn L = qfn_mul(&d, &inv);
    qfn_free(&d); qfn_free(&inv);
    return L;
}

// a/b as a rational constant, if it is one
static bool constant_ratio(const QFn* a, const QFn* b, Rational* r) {
    QFn ib = qfn_inv(b);
    QFn q = qfn_mul(a, &ib);
    bool c = polynomial_degree(q.num) <= 0 && polynomial_degree(q.den) == 0;
    if (c) poly_coeff(r, q.num, 0);
    qfn_free(&ib); qfn_free(&q);
    return c;
}

static QFn qfn_pow(const QFn* a, long k) {
    QFn base = k < 0 ? qfn_inv(a) : qfn_copy(a), r = { poly_one(), poly_one() };
    for (long i = 0; i < labs(k); i++) { QFn p = qfn_mul(&r, &base); qfn_free(&r); r = p; }
    qfn_free(&base);
    return r;
}

static bool smaller_abs(const Rational* a, const Rational* b) {
    Rational x, y;
    rat_init(&x); rat_init(&y);
    rat_abs(&x, a); rat_abs(&y, b);
    bool lt = rat_cmp(&x, &y) < 0;
    rat_clear(&x); rat_clear(&y);
    return lt;
}

// θ = ln(u) for the argument u with the smallest |ln u / ln w|; false unless every argument is
// exactly u^k for an integer k (ln(2x) = ln 2 + ln x brings in a new constant)
static bool make_log_tower(const Args* A, Tower* t) {
    QFn L0 = log_deriv(&A->w[0]);
    Rational* r = malloc((size_t)A->nw * sizeof(Rational));
    int base = 0;
    bool ok = true;
    for (int i = 0; i < A->nw; i++) {
        rat_init(&r[i]);
        QFn Li = log_deriv(&A->w[i]);
        ok = ok && constant_ratio(&Li, &L0, &r[i]);
        qfn_free(&Li);
        if (ok && smaller_abs(&r[i], &r[base])) base = i;
    }
    for (int i = 0; ok && i < A->nw; i++) {
        Rational k;
        rat_init(&k);
        rat_div(&k, &r[i], &r[base]);
        ok = mpz_cmp_ui(mpq_denref(k.val), 1) == 0 && mpz_cmpabs_ui(mpq_numref(k.val), MAX_EXP_POWER) <= 0;
        if (ok) {
            QFn p = qfn_pow(&A->w[base], mpz_get_si(mpq_numref(k.val)));
            QFn d = qfn_add(&p, &A->w[i], -1);
            ok = qfn_is_zero(&d);
            qfn_free(&p); qfn_free(&d);
        }
        rat_clear(&k);
    }
    if (ok) {                                             // orient as ln(x), not ln(1/x)
        QFn u = polynomial_degree(A->w[base].num) < polynomial_degree(A->w[base].den) ? qfn_inv(&A->w[base]) : qfn_copy(&A->w[base]);
        QFn L = log_deriv(&u);
        *t = (Tower){ TOWER_LOG, NULL, u.num, u.den, L.num, L.den, false };
    }
    for (int i = 0; i < A->nw; i++) rat_clear(&r[i]);
    free(r);
    qfn_free(&L0);
    return ok;
}

// the k with ln(w) = k·θ, for an argument w of the tower (an integer by construction)
static long log_multiple(const Tower* t, const QFn* w) {
    QFn Lw = log_deriv(w), L = tower_L(t);
    Rational k;
    rat_init(&k);
    constant_ratio(&Lw, &L, &k);
    long m = mpz_get_si(mpq_numref(k.val));
    rat_clear(&k);
    qfn_free(&Lw); qfn_free(&L);
    return m;
}

// --- expression -> fraction in θ ---

typedef struct {
    const Tower* t;
    const char* unsupported;
} Conv;

static KFrac conv(const Expr* e, Conv* C);

static KFrac kf_empty(void) { return (KFrac){ lp_empty(), lp_empty() }; }
static KFrac kf_poly(LPoly n) { return (KFrac){ n, kp_one() }; }  // takes n
static bool kf_is_poly(const KFrac* f) { return kp_is_one(&f->d); }

// 1/f for f ≠ 0
static KFrac kf_inv(const KFrac* f, const Conv* C) { return kfrac_inv(f, C->t->kind); }

static KFrac kf_mul(const KFrac* a, const KFrac* b, const Conv* C) {
    if (kf_is_poly(a) && kf_is_poly(b)) return kf_poly(lpoly_mul(&a->n, &b->n));
    return kfrac_mul(a, b, C->t->kind);
}

static KFrac kf_add(const KFrac* a, const KFrac* b, Conv* C) {
    if (kf_is_poly(a) && kf_is_poly(b)) return kf_poly(lpoly_add(&a->n, &b->n, 1));
    return kfrac_add(a, b, C->t->kind);
}

static KFrac kf_neg(const KFrac* f) {
    LPoly z = lp_empty();
    return (KFrac){ lpoly_add(&z, &f->n, -1), kp_copy(&f->d) };
}

static KFrac conv_pow(const Expr* e, Conv* C) {
    Expr* n = canonicalize(e->binop.right);
    bool integer = n && n->type == EXPR_POLY && polynomial_degree(n->poly) <= 0;
    long k = 0;
    if (integer) {
        Rational c;
        rat_init(&c);
        poly_coeff(&c, n->poly, 0);
        integer = mpz_cmp_ui(mpq_denref(c.val), 1) == 0 && mpz_cmpabs_ui(mpq_numref(c.val), MAX_EXP_POWER) <= 0;
        if (integer) k = mpz_get_si(mpq_numref(c.val));
        rat_clear(&c);
    }
    free_expr(n);
    if (!integer) { C->unsupported = "exponents must be small integers"; return kf_empty(); }
    KFrac b = conv(e->binop.left, C);
    if (C->unsupported) { kfrac_free(&b); return kf_empty(); }
    if (k < 0 && b.n.n == 0) { kfrac_free(&b); C->unsupported = "division by zero"; return kf_empty(); }
    if (k < 0) { KFrac inv = kf_inv(&b, C); kfrac_free(&b); b = inv; k = -k; }
    KFrac r = kf_poly(kp_one());
    for (long i = 0; i < k; i++) { KFrac p = kf_mul(&r, &b, C); kfrac_free(&r); r = p; }
    kfrac_free(&b);
    return r;
}

static KFrac conv_binop(const Expr* e, Conv* C) {
    if (e->binop.op == OP_POW) return conv_pow(e, C);
    KFrac l = conv(e->binop.left, C), r = conv(e->binop.right, C), out = kf_empty();
    if (!C->unsupported) switch (e->binop.op) {
        case OP_ADD: out = kf_add(&l, &r, C); break;
        case OP_SUB: { KFrac m = kf_neg(&r); out = kf_add(&l, &m, C); kfrac_free(&m); break; }
        case OP_MUL: out = kf_mul(&l, &r, C); break;
        default:
            if (r.n.n == 0) C->unsupported = "division by zero";
            else { KFrac inv = kf_inv(&r, C); out = kf_mul(&l, &inv, C); kfrac_free(&inv); }
    }
    kfrac_free(&l); kfrac_free(&r);
    return out;
}

// ln(u^k) = k·θ (collect checked u; ln(1) = 0)
static LPoly conv_ln(const Expr* e, const Conv* C) {
    Expr* c = canonicalize(e->ln_arg);
    QFn w = c->type == EXPR_POLY ? qfn_poly(c->poly) : qfn_from(poly_copy(c->ratfn.num), poly_copy(c->ratfn.den));
    free_expr(c);
    bool constant = polynomial_degree(w.num) <= 0 && polynomial_degree(w.den) == 0;
    Rational k = rat_from_int(constant ? 0 : log_multiple(C->t, &w));
    qfn_free(&w);
    QFn kq = { poly_constant(&k), poly_one() };
    rat_clear(&k);
    return lpoly_term(kq, 1);
}

static KFrac conv(const Expr* e, Conv* C) {
    if (C->unsupported) return kf_empty();
    if (!expr_has_exp(e) && !expr_has_ln(e)) {           // a rational function of x
        Expr* c = canonicalize((Expr*)e);                 // canonicalize does not modify e
        if (!c) { C->unsupported = "not a rational function of x and exp"; return kf_empty(); }
        QFn a = c->type == EXPR_POLY ? qfn_poly(c->poly) : qfn_from(poly_copy(c->ratfn.num), poly_copy(c->ratfn.den));
        free_expr(c);
        return kf_poly(lpoly_term(a, 0));
    }
    switch (e->type) {
        case EXPR_BINOP: return conv_binop(e, C);
        case EXPR_NEG: { KFrac f = conv(e->operand, C), r = kf_neg(&f); kfrac_free(&f); return r; }
        case EXPR_EXP: {
            Expr* c = canonicalize(e->fn_arg);           // a polynomial (collect checked)
            int k = polynomial_degree(c->poly) > 0 ? power_of(C->t, c->poly) : 0;
            free_expr(c);
            return kf_poly(lpoly_term((QFn){ poly_one(), poly_one() }, k));
        }
        case EXPR_LN: return kf_poly(conv_ln(e, C));
        default:
            C->unsupported = "only +, -, *, /, integer powers, exp and ln are supported";
            return kf_empty();
    }
}

ExpForm exp_form(const Expr* e) {
    ExpForm F = { { TOWER_EXP, NULL, NULL, NULL, NULL, NULL, false }, lp_empty(), kp_one(), NULL };
    Args A = { NULL, 0, NULL, 0, NULL };
    collect(e, &A);
    F.unsupported = A.unsupported;
    if (!F.unsupported && A.n > 0 && A.nw > 0) F.unsupported = "exp and ln together (two towers) are not supported yet";
    if (!F.unsupported && A.nw > 0 && !make_log_tower(&A, &F.t))
        F.unsupported = "ln of arguments that are not integer powers of one function (such as ln(x)·ln(x+1) or ln(2x)) is not supported yet";
    if (!F.unsupported && A.n > 0 && !make_tower(&A, &F.t))
        F.unsupported = "a denominator with exponentials of independent arguments (such as exp(x)/(exp(x) + exp(x^2))) is not supported yet";
    if (!F.unsupported && F.t.kind == TOWER_EXP && !F.t.v) F.t.v = poly_one(); // no θ: any v works
    if (!F.unsupported) {
        Conv C = { &F.t, NULL };
        KFrac q = conv(e, &C);
        F.unsupported = C.unsupported;
        if (!F.unsupported) { lpoly_free(&F.den); F.f = q.n; F.den = q.d; }
        else kfrac_free(&q);
    }
    for (int i = 0; i < A.n; i++) free_polynomial(A.u[i]);
    for (int i = 0; i < A.nw; i++) qfn_free(&A.w[i]);
    free(A.u); free(A.w);
    return F;
}

void exp_form_free(ExpForm* F) {
    lpoly_free(&F->f);
    lpoly_free(&F->den);
    Polynomial* parts[] = { F->t.v, F->t.un, F->t.ud, F->t.Ln, F->t.Ld };
    for (int i = 0; i < 5; i++) free_polynomial(parts[i]);
    F->t = (Tower){ TOWER_EXP, NULL, NULL, NULL, NULL, NULL, false };
}
