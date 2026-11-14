// Q(x)(θ), θ = exp(v): rational functions, Laurent polynomials in θ, and the conversion
// from an expression tree (which also finds v).
#include "expfield.h"
#include "modgcd.h"
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

static QFn qfn_poly(const Polynomial* p) { return (QFn){ poly_copy(p), poly_one() }; }

// a ± b
static QFn qfn_add(const QFn* a, const QFn* b, int sign) {
    Polynomial *l = poly_multiply(a->num, b->den), *r = poly_multiply(b->num, a->den);
    Polynomial* n = sign > 0 ? poly_add(l, r) : poly_subtract(l, r);
    free_polynomial(l); free_polynomial(r);
    return qfn_from(n, poly_multiply(a->den, b->den));
}

static QFn qfn_mul(const QFn* a, const QFn* b) {
    return qfn_from(poly_multiply(a->num, b->num), poly_multiply(a->den, b->den));
}

// 1/a for a ≠ 0
static QFn qfn_inv(const QFn* a) { return qfn_from(poly_copy(a->den), poly_copy(a->num)); }

// (n/d)' = (n'd - nd')/d²
static QFn qfn_derivative(const QFn* a) {
    Polynomial *dn = poly_derivative(a->num), *dd = poly_derivative(a->den);
    Polynomial *l = poly_multiply(dn, a->den), *r = poly_multiply(a->num, dd);
    Polynomial* n = poly_subtract(l, r);
    free_polynomial(dn); free_polynomial(dd); free_polynomial(l); free_polynomial(r);
    return qfn_from(n, poly_multiply(a->den, a->den));
}

Expr* qfn_expr(const QFn* a) {
    if (polynomial_degree(a->den) == 0) return create_poly_expr(poly_copy(a->num));
    return create_rational_fn(poly_copy(a->num), poly_copy(a->den));
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
static LPoly lp_term(QFn a, int k) {
    LPoly f = { k, 1, malloc(sizeof(QFn)) };
    f.a[0] = a;
    lp_trim(&f);
    return f;
}

static LPoly lp_add(const LPoly* f, const LPoly* g, int sign) {
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

static LPoly lp_mul(const LPoly* f, const LPoly* g) {
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

LPoly lpoly_derivative(const LPoly* f, const Tower* t) {
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
    LPoly d = lp_add(f, g, -1);
    bool z = d.n == 0;
    lpoly_free(&d);
    return z;
}

// exp(k·v)
static Expr* theta_power(const Tower* t, int k) {
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
    if (minus_one) e = create_neg(e);
    else if (!(constant && rat_is_one(&c))) e = create_mul(create_poly_expr(poly_copy(a->num)), e);
    rat_clear(&c);
    if (polynomial_degree(a->den) > 0) e = create_div(e, create_poly_expr(poly_copy(a->den)));
    return e;
}

Expr* lpoly_expr(const LPoly* f, const Tower* t) {
    Expr* acc = NULL;
    for (int i = 0; i < f->n; i++) {
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
    const char* unsupported;
} Args;

static void note(Args* A, const char* why) { if (!A->unsupported) A->unsupported = why; }

static void collect(const Expr* e, Args* A) {
    switch (e->type) {
        case EXPR_BINOP: collect(e->binop.left, A); collect(e->binop.right, A); return;
        case EXPR_NEG: collect(e->operand, A); return;
        case EXPR_LN: note(A, "ln is not supported yet"); return;
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

// --- expression -> Laurent polynomial ---

typedef struct {
    const Tower* t;
    const char* unsupported;
} Conv;

static LPoly conv(const Expr* e, Conv* C);

// a single term a·θ^k (the only divisors allowed)
static bool is_monomial(const LPoly* f) { return f->n == 1; }

static LPoly conv_pow(const Expr* e, Conv* C) {
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
    if (!integer) { C->unsupported = "exponents must be small integers"; return lp_empty(); }
    LPoly b = conv(e->binop.left, C);
    if (C->unsupported) { lpoly_free(&b); return lp_empty(); }
    if (k < 0 && !is_monomial(&b)) { lpoly_free(&b); C->unsupported = "division by a sum involving exp (such as 1/(1 + exp(x))) is not supported yet"; return lp_empty(); }
    if (k < 0) {                                          // (a·θ^j)^-1 = (1/a)·θ^-j
        LPoly inv = lp_term(qfn_inv(&b.a[0]), -b.lo);
        lpoly_free(&b);
        b = inv;
        k = -k;
    }
    LPoly r = lp_term((QFn){ poly_one(), poly_one() }, 0);
    for (long i = 0; i < k; i++) { LPoly p = lp_mul(&r, &b); lpoly_free(&r); r = p; }
    lpoly_free(&b);
    return r;
}

static LPoly conv_binop(const Expr* e, Conv* C) {
    if (e->binop.op == OP_POW) return conv_pow(e, C);
    LPoly l = conv(e->binop.left, C), r = conv(e->binop.right, C), out = lp_empty();
    if (!C->unsupported) switch (e->binop.op) {
        case OP_ADD: out = lp_add(&l, &r, 1); break;
        case OP_SUB: out = lp_add(&l, &r, -1); break;
        case OP_MUL: out = lp_mul(&l, &r); break;
        default:
            if (r.n == 0) C->unsupported = "division by zero";
            else if (!is_monomial(&r)) C->unsupported = "division by a sum involving exp (such as 1/(1 + exp(x))) is not supported yet";
            else { LPoly inv = lp_term(qfn_inv(&r.a[0]), -r.lo); out = lp_mul(&l, &inv); lpoly_free(&inv); }
    }
    lpoly_free(&l); lpoly_free(&r);
    return out;
}

static LPoly conv(const Expr* e, Conv* C) {
    if (C->unsupported) return lp_empty();
    if (!expr_has_exp(e)) {                               // a rational function of x
        Expr* c = canonicalize((Expr*)e);                 // canonicalize does not modify e
        if (!c) { C->unsupported = "not a rational function of x and exp"; return lp_empty(); }
        QFn a = c->type == EXPR_POLY ? qfn_poly(c->poly) : qfn_from(poly_copy(c->ratfn.num), poly_copy(c->ratfn.den));
        free_expr(c);
        return lp_term(a, 0);
    }
    switch (e->type) {
        case EXPR_BINOP: return conv_binop(e, C);
        case EXPR_NEG: {
            LPoly f = conv(e->operand, C), z = lp_empty();
            LPoly r = lp_add(&z, &f, -1);
            lpoly_free(&f);
            return r;
        }
        case EXPR_EXP: {
            Expr* c = canonicalize(e->fn_arg);           // a polynomial (collect checked)
            int k = polynomial_degree(c->poly) > 0 ? power_of(C->t, c->poly) : 0;
            free_expr(c);
            return lp_term((QFn){ poly_one(), poly_one() }, k);
        }
        default:
            C->unsupported = "only +, -, *, /, integer powers and exp are supported";
            return lp_empty();
    }
}

ExpForm exp_form(const Expr* e) {
    ExpForm F = { { NULL }, lp_empty(), NULL };
    Args A = { NULL, 0, NULL };
    collect(e, &A);
    F.unsupported = A.unsupported;
    if (!F.unsupported && A.n > 0 && !make_tower(&A, &F.t))
        F.unsupported = "exponentials of independent arguments (such as exp(x)·exp(x^2)) are not supported yet";
    if (!F.unsupported && !F.t.v) F.t.v = poly_one();     // no θ: exp(0) only; any v works
    if (!F.unsupported) {
        Conv C = { &F.t, NULL };
        F.f = conv(e, &C);
        F.unsupported = C.unsupported;
    }
    for (int i = 0; i < A.n; i++) free_polynomial(A.u[i]);
    free(A.u);
    return F;
}

void exp_form_free(ExpForm* F) {
    lpoly_free(&F->f);
    free_polynomial(F->t.v);
    F->t.v = NULL;
}
