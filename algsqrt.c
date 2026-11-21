// Q(x)(y), y = sqrt(U): conversion from expressions, arithmetic, derivative, output, and the
// integral of a + b·y through a rational parametrization of y² = U (deg U <= 2).
#include "algsqrt.h"
#include "expfield.h"
#include "kpoly.h"
#include "factor.h"
#include "residue.h"
#include "sqrtlog.h"
#include "partial_frac.h"
#include <math.h>
#include <stdlib.h>

// a + b·y
typedef struct {
    QFn a, b;
} SqFn;

// the radicand U (NULL until the first square root), and why a conversion failed
typedef struct {
    Polynomial* U;
    const char* why;
} SqCtx;

// --- arithmetic ---

static Polynomial* poly_x(void) {
    Polynomial* p = create_polynomial();
    Rational one = rat_one_val();
    add_term(p, &one, 1);
    rat_clear(&one);
    return p;
}

// c2·t² + c1·t + c0 (any of them 0)
static Polynomial* quad_t(const Rational* c2, const Rational* c1, const Rational* c0) {
    Polynomial* p = create_polynomial();
    add_term(p, c2, 2); add_term(p, c1, 1); add_term(p, c0, 0);
    return p;
}

static SqFn sq_zero(void) { return (SqFn){ qfn_zero(), qfn_zero() }; }
static SqFn sq_rat(QFn a) { return (SqFn){ a, qfn_zero() }; }      // takes a
static SqFn sq_copy(const SqFn* f) { return (SqFn){ qfn_copy(&f->a), qfn_copy(&f->b) }; }
static void sq_free(SqFn* f) { qfn_free(&f->a); qfn_free(&f->b); }
static bool sq_is_zero(const SqFn* f) { return qfn_is_zero(&f->a) && qfn_is_zero(&f->b); }

static SqFn sq_add(const SqFn* f, const SqFn* g, int sign) {
    return (SqFn){ qfn_add(&f->a, &g->a, sign), qfn_add(&f->b, &g->b, sign) };
}

// U as an element of Q(x), or 0 before any root (then every b is 0)
static QFn u_of(const SqCtx* c) { return c->U ? qfn_poly(c->U) : qfn_zero(); }

// (a1 + b1·y)(a2 + b2·y) = (a1·a2 + b1·b2·U) + (a1·b2 + a2·b1)·y
static SqFn sq_mul(const SqCtx* c, const SqFn* f, const SqFn* g) {
    QFn u = u_of(c), aa = qfn_mul(&f->a, &g->a), bb = qfn_mul(&f->b, &g->b), ubb = qfn_mul(&u, &bb);
    QFn ab = qfn_mul(&f->a, &g->b), ba = qfn_mul(&f->b, &g->a);
    SqFn r = { qfn_add(&aa, &ubb, 1), qfn_add(&ab, &ba, 1) };
    qfn_free(&u); qfn_free(&aa); qfn_free(&bb); qfn_free(&ubb); qfn_free(&ab); qfn_free(&ba);
    return r;
}

// 1/(a + b·y) = (a - b·y)/(a² - b²·U) for f ≠ 0 (the norm vanishes only at 0: U is not a square)
static SqFn sq_inv(const SqCtx* c, const SqFn* f) {
    QFn u = u_of(c), aa = qfn_mul(&f->a, &f->a), bb = qfn_mul(&f->b, &f->b), ubb = qfn_mul(&u, &bb);
    QFn n = qfn_add(&aa, &ubb, -1), in = qfn_inv(&n), z = qfn_zero(), nb = qfn_add(&z, &f->b, -1);
    SqFn r = { qfn_mul(&f->a, &in), qfn_mul(&nb, &in) };
    qfn_free(&u); qfn_free(&aa); qfn_free(&bb); qfn_free(&ubb); qfn_free(&n); qfn_free(&in); qfn_free(&z); qfn_free(&nb);
    return r;
}

// f^k for any integer k (f ≠ 0 when k < 0), by squaring
static SqFn sq_pow(const SqCtx* c, const SqFn* f, long k) {
    SqFn r = sq_rat(qfn_const(1, 1)), base = k < 0 ? sq_inv(c, f) : sq_copy(f);
    for (unsigned long e = (unsigned long)(k < 0 ? -k : k); e; e >>= 1) {
        if (e & 1) { SqFn t = sq_mul(c, &r, &base); sq_free(&r); r = t; }
        if (e > 1) { SqFn t = sq_mul(c, &base, &base); sq_free(&base); base = t; }
    }
    sq_free(&base);
    return r;
}

// D(a + b·y) = a' + (b' + b·U'/(2U))·y
static SqFn sq_deriv(const SqCtx* c, const SqFn* f) {
    SqFn r = { qfn_derivative(&f->a), qfn_derivative(&f->b) };
    if (!c->U || qfn_is_zero(&f->b)) return r;
    Rational two = rat_from_int(2);
    QFn half = qfn_from(poly_derivative(c->U), poly_scale(c->U, &two)), bh = qfn_mul(&f->b, &half), s = qfn_add(&r.b, &bh, 1);
    qfn_free(&r.b); qfn_free(&half); qfn_free(&bh);
    r.b = s;
    rat_clear(&two);
    return r;
}

// --- the radicand: sqrt(g) = r·sqrt(U), U = m·w unique up to positive squares ---

static bool rational_square(const Rational* q) {
    return mpq_sgn(q->val) >= 0 && mpz_perfect_square_p(mpq_numref(q->val)) && mpz_perfect_square_p(mpq_denref(q->val));
}

// every factor (and the content) to an even power
static bool poly_is_square(const Polynomial* p) {
    Factorization* f = factorize(p);
    bool sq = rational_square(&f->content);
    for (int i = 0; sq && i < f->count; i++) sq = f->multiplicities[i] % 2 == 0;
    free_factorization(f);
    return sq;
}

static bool same_poly(const Polynomial* p, const Polynomial* q) {
    Polynomial* d = poly_subtract(p, q);
    bool z = poly_is_zero(d);
    free_polynomial(d);
    return z;
}

// ±m·w for g = c·w (w primitive with lc(w) > 0, |c| = r²·m), the sign that of lc(g), and r
static Polynomial* normal_radicand(const Polynomial* g, Rational* r) {
    Rational cont = poly_content(g), lc, lw;
    rat_init(&lc); rat_init(&lw);
    poly_coeff(&lc, g, polynomial_degree(g));
    Polynomial* w = poly_primitive_part(g);
    poly_coeff(&lw, w, polynomial_degree(w));
    if (rat_is_negative(&lw)) { Rational m1 = rat_from_int(-1); Polynomial* n = poly_scale(w, &m1); free_polynomial(w); w = n; rat_clear(&m1); }
    rat_abs(&cont, &cont);
    mpz_t m;
    mpz_init(m);
    rat_sqrt_split(&cont, r, m);
    Rational k;
    rat_init(&k);
    mpq_set_z(k.val, m);
    if (rat_is_negative(&lc)) rat_neg(&k, &k);
    Polynomial* U = poly_scale(w, &k);
    free_polynomial(w);
    rat_clear(&cont); rat_clear(&lc); rat_clear(&lw); rat_clear(&k);
    mpz_clear(m);
    return U;
}

// sqrt(g) for a polynomial g: r·y with U fixed by the first root, or a rational r when g = r²
static bool radical(SqCtx* c, const Polynomial* g, SqFn* out) {
    if (poly_is_zero(g)) { *out = sq_zero(); return true; }
    Rational r, u0;
    rat_init(&r); rat_init(&u0);
    Polynomial* U = normal_radicand(g, &r);
    poly_coeff(&u0, U, 0);
    bool constant = polynomial_degree(U) == 0, ok = true;
    if (constant && rat_is_one(&u0)) *out = sq_rat(qfn_rat(&r));                  // sqrt(4/9) = 2/3
    else if (constant && rat_is_negative(&u0)) { ok = false; c->why = "the square root of a negative number is not real"; }
    else if (!constant && poly_is_square(U)) { ok = false; c->why = "the square root of a perfect square is an absolute value, which is not supported"; }
    else if (c->U && !same_poly(c->U, U)) { ok = false; c->why = "two different square roots in one expression are not supported yet"; }
    else {
        if (!c->U) c->U = poly_copy(U);
        *out = (SqFn){ qfn_zero(), qfn_rat(&r) };
    }
    free_polynomial(U);
    rat_clear(&r); rat_clear(&u0);
    return ok;
}

// --- expressions -> Q(x)(y) ---

// e as a rational constant (an exponent); false when it depends on x or is not a rational function
static bool constant_of(const Expr* e, Rational* k) {
    Expr* c = canonicalize((Expr*)e);
    bool ok = c && c->type == EXPR_POLY && polynomial_degree(c->poly) <= 0;
    if (ok) poly_coeff(k, c->poly, 0);
    free_expr(c);
    return ok;
}

static bool conv(const Expr* e, SqCtx* c, SqFn* out);

// a polynomial operand (the argument of a square root)
static bool conv_poly(const Expr* e, SqCtx* c, Polynomial** p) {
    SqFn f;
    if (!conv(e, c, &f)) return false;
    bool ok = qfn_is_zero(&f.b) && polynomial_degree(f.a.den) == 0;
    if (ok) *p = poly_copy(f.a.num);                   // the denominator is monic: 1
    else c->why = "a square root of a rational function or of another root is not supported yet";
    sq_free(&f);
    return ok;
}

// base^k for k an integer or a half-integer n/2 (= sqrt(base)^n)
static bool conv_pow(const Expr* e, SqCtx* c, SqFn* out) {
    Rational k;
    rat_init(&k);
    bool ok = constant_of(e->binop.right, &k);
    long num = ok ? mpz_get_si(mpq_numref(k.val)) : 0, den = ok ? mpz_get_si(mpq_denref(k.val)) : 0;
    ok = ok && (den == 1 || den == 2) && mpz_fits_slong_p(mpq_numref(k.val)) && labs(num) <= 256;
    if (!ok) c->why = "exponents must be integers or halves of integers (up to 256)";
    SqFn base;
    Polynomial* g = NULL;
    if (ok && den == 2) ok = conv_poly(e->binop.left, c, &g) && radical(c, g, &base);
    else if (ok) ok = conv(e->binop.left, c, &base);
    if (ok && num < 0 && sq_is_zero(&base)) { ok = false; c->why = "division by zero"; sq_free(&base); }
    if (ok) { *out = sq_pow(c, &base, num); sq_free(&base); }
    free_polynomial(g);
    rat_clear(&k);
    return ok;
}

static bool conv_binop(const Expr* e, SqCtx* c, SqFn* out) {
    if (e->binop.op == OP_POW) return conv_pow(e, c, out);
    SqFn l, r;
    if (!conv(e->binop.left, c, &l)) return false;
    if (!conv(e->binop.right, c, &r)) { sq_free(&l); return false; }
    bool ok = true;
    switch (e->binop.op) {
        case OP_ADD: *out = sq_add(&l, &r, 1); break;
        case OP_SUB: *out = sq_add(&l, &r, -1); break;
        case OP_MUL: *out = sq_mul(c, &l, &r); break;
        default:                                          // OP_DIV
            if (sq_is_zero(&r)) { ok = false; c->why = "division by zero"; break; }
            SqFn ir = sq_inv(c, &r);
            *out = sq_mul(c, &l, &ir);
            sq_free(&ir);
    }
    sq_free(&l); sq_free(&r);
    return ok;
}

static bool conv(const Expr* e, SqCtx* c, SqFn* out) {
    switch (e->type) {
        case EXPR_CONST: *out = sq_rat(qfn_rat(&e->constant)); return true;
        case EXPR_VAR: { Polynomial* x = poly_x(); *out = sq_rat(qfn_poly(x)); free_polynomial(x); return true; }
        case EXPR_POLY: *out = sq_rat(qfn_poly(e->poly)); return true;
        case EXPR_RATIONAL_FN: *out = sq_rat(qfn_from(poly_copy(e->ratfn.num), poly_copy(e->ratfn.den))); return true;
        case EXPR_NEG: {
            SqFn f;
            if (!conv(e->operand, c, &f)) return false;
            SqFn z = sq_zero();
            *out = sq_add(&z, &f, -1);
            sq_free(&z); sq_free(&f);
            return true;
        }
        case EXPR_BINOP: return conv_binop(e, c, out);
        case EXPR_SQRT: {
            Polynomial* g = NULL;
            bool ok = conv_poly(e->fn_arg, c, &g) && radical(c, g, out);
            free_polynomial(g);
            return ok;
        }
        default:
            c->why = "a square root together with exp, ln, sin or cos is not supported yet";
            return false;
    }
}

// --- Q(x)(y) -> expressions ---

// b·sqrt(U): c·sqrt(U), p·sqrt(U), or (n·sqrt(U))/d
static Expr* times_root(const QFn* b, const Polynomial* U) {
    Rational u;
    rat_init(&u);
    poly_coeff(&u, U, 0);
    Expr* s = create_sqrt(polynomial_degree(U) == 0 ? create_const_rat(&u) : create_poly_expr(poly_copy(U))); // sqrt(2) a constant
    rat_clear(&u);
    Rational k;
    rat_init(&k);
    Expr* n;
    if (polynomial_degree(b->num) == 0) {                 // a constant numerator folds in front
        poly_coeff(&k, b->num, 0);
        n = rat_is_one(&k) ? s : create_mul(create_const_rat(&k), s);
    } else n = create_mul(create_poly_expr(poly_copy(b->num)), s);
    rat_clear(&k);
    return polynomial_degree(b->den) == 0 ? n : create_div(n, create_poly_expr(poly_copy(b->den)));
}

// a + b·y, over one denominator when a and b have a nonconstant one: (A + B·y)/D
static Expr* sq_expr(const SqFn* f, const Polynomial* U) {
    if (qfn_is_zero(&f->b)) return qfn_expr(&f->a);
    if (qfn_is_zero(&f->a)) return times_root(&f->b, U);
    if (polynomial_degree(f->a.den) == 0 && polynomial_degree(f->b.den) == 0) return create_add(qfn_expr(&f->a), times_root(&f->b, U));
    Polynomial* D = poly_lcm(f->a.den, f->b.den);              // both monic, so D is too
    QFn d = qfn_poly(D), A = qfn_mul(&f->a, &d), B = qfn_mul(&f->b, &d);   // polynomials now
    Expr* e = create_div(create_add(qfn_expr(&A), times_root(&B, U)), create_poly_expr(D));
    qfn_free(&d); qfn_free(&A); qfn_free(&B);
    return e;
}

// e converted, or false with c->why (U is freed by the caller)
static bool form_of(const Expr* e, SqCtx* c, SqFn* f) {
    *c = (SqCtx){ NULL, NULL };
    return conv(e, c, f);
}

Expr* simplify_sqrt(const Expr* e, const char** why) {
    SqCtx c;
    SqFn f;
    Expr* r = NULL;
    if (form_of(e, &c, &f)) { r = sq_expr(&f, c.U); sq_free(&f); }
    *why = c.why;
    free_polynomial(c.U);
    return r;
}

Expr* differentiate_sqrt(const Expr* e, const char** why) {
    SqCtx c;
    SqFn f;
    Expr* r = NULL;
    if (form_of(e, &c, &f)) {
        SqFn d = sq_deriv(&c, &f);
        r = sq_expr(&d, c.U);
        sq_free(&d); sq_free(&f);
    }
    *why = c.why;
    free_polynomial(c.U);
    return r;
}

// a sqrt node, or a power whose exponent is not an integer
bool expr_has_sqrt(const Expr* e) {
    switch (e->type) {
        case EXPR_SQRT: return true;
        case EXPR_NEG: return expr_has_sqrt(e->operand);
        case EXPR_LN: return expr_has_sqrt(e->ln_arg);
        case EXPR_ATAN: case EXPR_EXP: case EXPR_SIN: case EXPR_COS: return expr_has_sqrt(e->fn_arg);
        case EXPR_BINOP: {
            if (expr_has_sqrt(e->binop.left) || expr_has_sqrt(e->binop.right)) return true;
            if (e->binop.op != OP_POW) return false;
            Rational k;
            rat_init(&k);
            bool half = constant_of(e->binop.right, &k) && mpz_cmp_ui(mpq_denref(k.val), 2) == 0;
            rat_clear(&k);
            return half;
        }
        default: return false;
    }
}

// --- a rational parametrization x = X(t), y = Y(t) of y² = U, and its inverse t = T(x, y) ---

typedef struct {
    QFn X, Y;   // rational functions of t (written in the variable x)
    SqFn T;     // t as an element of Q(x)(y)
} Param;

static void param_free(Param* P) { qfn_free(&P->X); qfn_free(&P->Y); sq_free(&P->T); }

// X = Xn/D, Y = Yn/D (takes all three)
static void param_set(Param* P, Polynomial* Xn, Polynomial* Yn, Polynomial* D) {
    P->X = qfn_from(Xn, poly_copy(D));
    P->Y = qfn_from(Yn, D);
}

// U = b·x + c: t = y, X = (t² - c)/b, Y = b·t/b
static void param_linear(const Rational* q, Param* P) {
    Rational one = rat_one_val(), zero = rat_zero_val(), mc;
    rat_init(&mc);
    rat_neg(&mc, &q[0]);
    param_set(P, quad_t(&one, &zero, &mc), quad_t(&zero, &q[1], &zero), poly_constant(&q[1]));
    P->T = (SqFn){ qfn_zero(), qfn_const(1, 1) };
    rat_clear(&one); rat_clear(&zero); rat_clear(&mc);
}

// a = s²: y = t - s·x (Euler), X = (t² - c)/(2s·t + b), Y = (s·t² + b·t + s·c)/(2s·t + b), T = s·x + y
static void param_euler(const Rational* q, const Rational* s, Param* P) {
    Rational one = rat_one_val(), zero = rat_zero_val(), mc, s2, sc;
    rat_init(&mc); rat_init(&s2); rat_init(&sc);
    rat_neg(&mc, &q[0]);
    rat_add(&s2, s, s);
    rat_mul(&sc, s, &q[0]);
    param_set(P, quad_t(&one, &zero, &mc), quad_t(s, &q[1], &sc), quad_t(&zero, &s2, &q[1]));
    Polynomial* sx = quad_t(&zero, s, &zero);
    P->T = (SqFn){ qfn_poly(sx), qfn_const(1, 1) };
    free_polynomial(sx);
    rat_clear(&one); rat_clear(&zero); rat_clear(&mc); rat_clear(&s2); rat_clear(&sc);
}

// a point (x0, y0), y0² = U(x0): t = (y - y0)/(x - x0), with u1 = U'(x0)
//   X = (x0·t² - 2y0·t + u1 - a·x0)/(t² - a),  Y = (-y0·t² + u1·t - a·y0)/(t² - a)
static void param_point(const Rational* q, const Rational* pt, Param* P) {
    const Rational *x0 = &pt[0], *y0 = &pt[1], *a = &q[2];
    Rational one = rat_one_val(), zero = rat_zero_val(), u1, t, m2y, k, ma, my, may;
    rat_init(&u1); rat_init(&t); rat_init(&m2y); rat_init(&k); rat_init(&ma); rat_init(&my); rat_init(&may);
    rat_mul(&u1, a, x0); rat_add(&u1, &u1, &u1); rat_add(&u1, &u1, &q[1]);   // 2a·x0 + b
    rat_add(&m2y, y0, y0); rat_neg(&m2y, &m2y);
    rat_mul(&t, a, x0); rat_sub(&k, &u1, &t);                                // u1 - a·x0
    rat_neg(&ma, a); rat_neg(&my, y0); rat_mul(&may, &ma, y0);
    param_set(P, quad_t(x0, &m2y, &k), quad_t(&my, &u1, &may), quad_t(&one, &zero, &ma));
    Rational mx0;
    rat_init(&mx0);
    rat_neg(&mx0, x0);
    Polynomial* lin = quad_t(&zero, &one, &mx0);                             // x - x0
    QFn d = qfn_poly(lin), id = qfn_inv(&d), ny = qfn_rat(&my);
    P->T = (SqFn){ qfn_mul(&ny, &id), id };
    qfn_free(&d); qfn_free(&ny);
    free_polynomial(lin);
    rat_clear(&one); rat_clear(&zero); rat_clear(&u1); rat_clear(&t); rat_clear(&m2y); rat_clear(&k);
    rat_clear(&ma); rat_clear(&my); rat_clear(&may); rat_clear(&mx0);
}

#define POINT_HEIGHT 40 // rational points x0 = p/q with |p| + q up to this are tried

// the square root of a rational square
static void rat_sqrt_exact(Rational* r, const Rational* q) {
    mpz_sqrt(mpq_numref(r->val), mpq_numref(q->val));
    mpz_sqrt(mpq_denref(r->val), mpq_denref(q->val));
}

// a rational point (x0, y0 >= 0) on y² = U, small x0 first (x0 = 0 is tried first)
static bool rational_point(const Polynomial* U, Rational* pt) {
    Rational v;
    rat_init(&v);
    bool found = false;
    for (int h = 1; !found && h <= POINT_HEIGHT; h++)
        for (int qd = 1; !found && qd <= h; qd++)
            for (int sg = 1; !found && sg >= -1; sg -= 2) {
                int p = sg * (h - qd);
                if (p == 0 && (qd != 1 || sg < 0)) continue;      // 0 once
                mpq_set_si(pt[0].val, p, (unsigned long)qd);
                mpq_canonicalize(pt[0].val);
                if (mpz_cmp_ui(mpq_denref(pt[0].val), (unsigned long)qd) != 0) continue;   // not in lowest terms: seen before
                poly_eval_rat(&v, U, &pt[0]);
                if (rational_square(&v)) { rat_sqrt_exact(&pt[1], &v); found = true; }
            }
    rat_clear(&v);
    return found;
}

// U's coefficients q[k] of x^k, k <= 2
static void quad_coeffs(const Polynomial* U, Rational* q) {
    for (int k = 0; k < 3; k++) { rat_init(&q[k]); poly_coeff(&q[k], U, k); }
}

// a parametrization of y² = U for deg U = 1 or 2; false when U is a conic without a rational point
static bool param_of(const Polynomial* U, Param* P) {
    Rational q[3], s, pt[2];
    quad_coeffs(U, q);
    rat_init(&s); rat_init(&pt[0]); rat_init(&pt[1]);
    bool ok = true;
    if (polynomial_degree(U) == 1) param_linear(q, P);
    else if (rational_square(&q[2])) { rat_sqrt_exact(&s, &q[2]); param_euler(q, &s, P); }
    else if (rational_point(U, pt)) param_point(q, pt, P);
    else ok = false;
    for (int k = 0; k < 3; k++) rat_clear(&q[k]);
    rat_clear(&s); rat_clear(&pt[0]); rat_clear(&pt[1]);
    return ok;
}

// --- composition ---

// p(X) for X ∈ Q(t), by Horner
static QFn poly_at_q(const Polynomial* p, const QFn* X) {
    QFn acc = qfn_zero();
    Rational c;
    rat_init(&c);
    for (int k = polynomial_degree(p); k >= 0; k--) {
        poly_coeff(&c, p, k);
        QFn ck = qfn_rat(&c), m = qfn_mul(&acc, X), s = qfn_add(&m, &ck, 1);
        qfn_free(&acc); qfn_free(&ck); qfn_free(&m);
        acc = s;
    }
    rat_clear(&c);
    return acc;
}

// f(X) for f, X ∈ Q(t)
static QFn qfn_at_q(const QFn* f, const QFn* X) {
    QFn n = poly_at_q(f->num, X), d = poly_at_q(f->den, X), id = qfn_inv(&d), r = qfn_mul(&n, &id);
    qfn_free(&n); qfn_free(&d); qfn_free(&id);
    return r;
}

// p(T) for T ∈ Q(x)(y)
static SqFn poly_at_sq(const SqCtx* c, const Polynomial* p, const SqFn* T) {
    SqFn acc = sq_zero();
    Rational k;
    rat_init(&k);
    for (int e = polynomial_degree(p); e >= 0; e--) {
        poly_coeff(&k, p, e);
        SqFn ck = sq_rat(qfn_rat(&k)), m = sq_mul(c, &acc, T), s = sq_add(&m, &ck, 1);
        sq_free(&acc); sq_free(&ck); sq_free(&m);
        acc = s;
    }
    rat_clear(&k);
    return acc;
}

// f(T) for f ∈ Q(t), T ∈ Q(x)(y)
static SqFn qfn_at_sq(const SqCtx* c, const QFn* f, const SqFn* T) {
    SqFn n = poly_at_sq(c, f->num, T), d = poly_at_sq(c, f->den, T), id = sq_inv(c, &d), r = sq_mul(c, &n, &id);
    sq_free(&n); sq_free(&d); sq_free(&id);
    return r;
}

static bool sq_equal(const SqFn* f, const SqFn* g) {
    SqFn d = sq_add(f, g, -1);
    bool z = sq_is_zero(&d);
    sq_free(&d);
    return z;
}

// X(T) = x and Y(T) = y in Q(x)(y): then F(T) has derivative f(x, y) whenever F'(t) = f(X, Y)·X'
static bool param_inverts(const SqCtx* c, const Param* P) {
    Polynomial* x = poly_x();
    SqFn X = qfn_at_sq(c, &P->X, &P->T), Y = qfn_at_sq(c, &P->Y, &P->T), ex = sq_rat(qfn_poly(x)), ey = { qfn_zero(), qfn_const(1, 1) };
    bool ok = sq_equal(&X, &ex) && sq_equal(&Y, &ey);
    sq_free(&X); sq_free(&Y); sq_free(&ex); sq_free(&ey);
    free_polynomial(x);
    return ok;
}

// --- back from t to x ---

typedef struct {
    const SqCtx* c;
    const SqFn* T;
} Back;

static Expr* back(const Expr* F, const Back* B);

// Σ_k s[k](α)·t^k of a rational RootSum as Σ_i α^i·c_i(x, y), the c_i ∈ Q(x)(y)
static Expr* back_rootsum(const Expr* F, const Back* B) {
    int n = polynomial_degree(F->rootsum.q);
    Expr** c = malloc((size_t)n * sizeof(Expr*));
    Rational k;
    rat_init(&k);
    for (int i = 0; i < n; i++) {
        Polynomial* p = create_polynomial();                  // Σ_k [α^i]s[k]·t^k
        for (int e = 0; e <= F->rootsum.s_deg; e++) { poly_coeff(&k, F->rootsum.s[e], i); add_term(p, &k, e); }
        SqFn v = poly_at_sq(B->c, p, B->T);
        c[i] = sq_is_zero(&v) ? NULL : sq_expr(&v, B->c->U);
        sq_free(&v);
        free_polynomial(p);
    }
    rat_clear(&k);
    return create_rootsum_terms(poly_copy(F->rootsum.q), c);
}

// F(t) with t = T(x, y): rational pieces through Q(x)(y), everything else node by node
static Expr* back(const Expr* F, const Back* B) {
    switch (F->type) {
        case EXPR_VAR: return sq_expr(B->T, B->c->U);
        case EXPR_POLY: { SqFn v = poly_at_sq(B->c, F->poly, B->T); Expr* e = sq_expr(&v, B->c->U); sq_free(&v); return e; }
        case EXPR_RATIONAL_FN: {
            QFn f = qfn_from(poly_copy(F->ratfn.num), poly_copy(F->ratfn.den));
            SqFn v = qfn_at_sq(B->c, &f, B->T);
            Expr* e = sq_expr(&v, B->c->U);
            sq_free(&v); qfn_free(&f);
            return e;
        }
        case EXPR_BINOP: {
            Expr *l = back(F->binop.left, B), *r = back(F->binop.right, B);
            switch (F->binop.op) {
                case OP_ADD: return create_add(l, r);
                case OP_SUB: return create_sub(l, r);
                case OP_MUL: return create_mul(l, r);
                case OP_DIV: return create_div(l, r);
                default: return create_pow(l, r);
            }
        }
        case EXPR_NEG: return create_neg(back(F->operand, B));
        case EXPR_LN: return create_ln(back(F->ln_arg, B));
        case EXPR_ATAN: return create_atan(back(F->fn_arg, B));
        case EXPR_ROOTSUM: return F->rootsum.s ? back_rootsum(F, B) : expr_copy(F);
        default: return expr_copy(F);                         // constants: sqrt(2), cos(pi/5), ...
    }
}

// F's summands: those rational in t (no ln, atan, RootSum or irrational constant) summed in
// Q(t) and substituted once, so they come back as one normal form a + b·y; the rest term by term
typedef struct {
    QFn rat;
    Expr* rest;
} Split;

static void split_terms(const Expr* F, int sign, Split* S, const Back* B) {
    if (F->type == EXPR_BINOP && (F->binop.op == OP_ADD || F->binop.op == OP_SUB)) {
        split_terms(F->binop.left, sign, S, B);
        split_terms(F->binop.right, F->binop.op == OP_ADD ? sign : -sign, S, B);
        return;
    }
    if (F->type == EXPR_NEG) { split_terms(F->operand, -sign, S, B); return; }
    Expr* c = canonicalize((Expr*)F);
    if (c && (c->type == EXPR_POLY || c->type == EXPR_RATIONAL_FN)) {
        QFn q = c->type == EXPR_POLY ? qfn_poly(c->poly) : qfn_from(poly_copy(c->ratfn.num), poly_copy(c->ratfn.den));
        QFn s = qfn_add(&S->rat, &q, sign);
        qfn_free(&S->rat); qfn_free(&q);
        S->rat = s;
    } else {
        Expr* t = back(F, B);
        S->rest = expr_append_sum(S->rest, sign > 0 ? t : create_neg(t));
    }
    free_expr(c);
}

// the whole of F(T): the rational part as one element of Q(x)(y), then the other terms
static Expr* back_sum(const Expr* F, const Back* B) {
    Split S = { qfn_zero(), NULL };
    split_terms(F, 1, &S, B);
    Expr* e = NULL;
    if (!qfn_is_zero(&S.rat)) {
        SqFn v = qfn_at_sq(B->c, &S.rat, B->T);
        Rational k;
        rat_init(&k);
        if (qfn_constant(&v.a, &k)) { qfn_free(&v.a); v.a = qfn_zero(); }   // an additive constant: dropped
        rat_clear(&k);
        if (!sq_is_zero(&v)) e = sq_expr(&v, B->c->U);
        sq_free(&v);
    }
    qfn_free(&S.rat);
    e = S.rest ? expr_append_sum(e, S.rest) : e;
    return e ? e : create_const(0);
}

// --- integration ---

// --- conics without a rational point: ∫ N/(U^m·y) by reduction (U = a·x² + b·x + c) ---

// b·y = N/(U^m·y) for b = P/U^k: (N, m) = (P·U, 0) when k = 0, else (P, k - 1); false when den(b)
// is not a power of U (the monic den is (U/lc)^k, so P = num·lc^k)
static bool as_reduction(const Polynomial* U, const QFn* b, Polynomial** N, int* m) {
    int k = polynomial_degree(b->den) / 2;
    Rational lc, lk;
    rat_init(&lc); rat_init(&lk);
    poly_coeff(&lc, U, 2);
    Polynomial* Um = poly_copy(U);
    poly_make_monic(Um);
    Polynomial* Uk = poly_pow(Um, (unsigned)k);
    bool ok = polynomial_degree(b->den) == 2 * k && same_poly(Uk, b->den);
    if (ok) {
        rat_pow_int(&lk, &lc, k);
        Polynomial* P = poly_scale(b->num, &lk);
        *N = k == 0 ? poly_multiply(P, U) : poly_copy(P);
        *m = k == 0 ? 0 : k - 1;
        free_polynomial(P);
    }
    free_polynomial(Um); free_polynomial(Uk);
    rat_clear(&lc); rat_clear(&lk);
    return ok;
}

// Q'·U + h·Q·U' for the derivative of Q·y/U^m (h = 1/2 - m): D(Q·y/U^m) = (Q'U + h·Q·U')/(U^m·y)
static Polynomial* deriv_numerator(const Polynomial* U, const Polynomial* Q, const Rational* h) {
    Polynomial *dQ = poly_derivative(Q), *dU = poly_derivative(U), *l = poly_multiply(dQ, U), *qu = poly_multiply(Q, dU);
    Polynomial *r = poly_scale(qu, h), *s = poly_add(l, r);
    free_polynomial(dQ); free_polynomial(dU); free_polynomial(l); free_polynomial(qu); free_polynomial(r);
    return s;
}

// one step for m >= 1: N/(U^m·y) = D(Q·y/U^m) + N1/(U^(m-1)·y), with Q ≡ N·U'^-1/h mod U so that
// U divides N - (Q'U + h·Q·U'); returns N1
static Polynomial* reduce_step(const Polynomial* U, const Polynomial* N, int m, Polynomial** Q) {
    Rational h;
    rat_init(&h);
    mpq_set_si(h.val, 1 - 2 * m, 2);
    Polynomial* dU = poly_derivative(U);
    PolyBezout bz = poly_xgcd(dU, U);                         // s·U' + t·U = 1 (U square-free)
    Polynomial *ns = poly_multiply(N, bz.s), *q = poly_mod(ns, U), *ih = create_polynomial();
    Rational inv;
    rat_init(&inv);
    mpq_inv(inv.val, h.val);
    *Q = poly_scale(q, &inv);
    Polynomial *D = deriv_numerator(U, *Q, &h), *rest = poly_subtract(N, D), *N1 = poly_divide(rest, U);
    free_polynomial(dU); free_polynomial(ns); free_polynomial(q); free_polynomial(ih); free_polynomial(D); free_polynomial(rest);
    free_bezout(&bz);
    rat_clear(&h); rat_clear(&inv);
    return N1;
}

// m = 0: N/y = D(Q·y) + λ/y with N = Q'U + Q·U'/2 + λ, Q solved from the top coefficient down
// (x^j in Q puts a·(j + 1) on x^(j+1)); returns Q
static Polynomial* reduce_base(const Polynomial* U, const Polynomial* N, Rational* lambda) {
    Rational a, c, half;
    rat_init(&a); rat_init(&c); rat_init(&half);
    poly_coeff(&a, U, 2);
    mpq_set_si(half.val, 1, 2);
    Polynomial *Q = create_polynomial(), *R = poly_copy(N);
    for (int j = polynomial_degree(N) - 1; j >= 0; j--) {
        poly_coeff(&c, R, j + 1);
        Rational k = rat_from_int(j + 1);
        rat_mul(&k, &k, &a);
        rat_div(&c, &c, &k);                                  // q_j
        Polynomial* t = create_polynomial();
        add_term(t, &c, j);
        Polynomial *D = deriv_numerator(U, t, &half), *R2 = poly_subtract(R, D), *Q2 = poly_add(Q, t);
        free_polynomial(t); free_polynomial(D); free_polynomial(R); free_polynomial(Q);
        R = R2; Q = Q2;
        rat_clear(&k);
    }
    poly_coeff(lambda, R, 0);                                 // R is the constant λ now
    free_polynomial(R);
    rat_clear(&a); rat_clear(&c); rat_clear(&half);
    return Q;
}

// λ·∫ 1/y for deg U = 2: (λ/√A)·atan(√A·(x + b/(2a))/y), A = -a, or (λ/√a)·ln|2a·x + b + 2√a·y|
static Expr* base_integral(const Polynomial* U, const Rational* lambda) {
    Rational a, b, A, r, c1, z0;
    rat_init(&a); rat_init(&b); rat_init(&A); rat_init(&r); rat_init(&c1); rat_init(&z0);
    poly_coeff(&a, U, 2); poly_coeff(&b, U, 1);
    bool neg = rat_is_negative(&a);
    rat_abs(&A, &a);
    mpz_t m;
    mpz_init(m);
    rat_sqrt_split(&A, &r, m);                                // sqrt|a| = r·sqrt(m)
    Rational mq;
    rat_init(&mq);
    mpq_set_z(mq.val, m);
    rat_mul(&c1, &r, &mq);
    mpq_inv(c1.val, c1.val);
    rat_mul(&c1, &c1, lambda);                                // λ/sqrt|a| = c1·sqrt(m)
    Expr* y = create_sqrt(create_poly_expr(poly_copy(U))), *body;
    Polynomial* lin = create_polynomial();
    if (neg) {                                                // atan(r·sqrt(m)·(x + b/(2a))/y)
        rat_add(&z0, &a, &a);
        rat_div(&z0, &b, &z0);
        Rational one = rat_one_val();
        add_term(lin, &one, 1); add_term(lin, &z0, 0);
        rat_clear(&one);
        body = create_atan(create_div(surd_times(&r, m, create_poly_expr(lin)), y));
    } else {                                                  // ln|2a·x + b + 2r·sqrt(m)·y|
        Rational a2, r2;
        rat_init(&a2); rat_init(&r2);
        rat_add(&a2, &a, &a); rat_add(&r2, &r, &r);
        add_term(lin, &a2, 1); add_term(lin, &b, 0);
        body = create_ln(create_add(create_poly_expr(lin), surd_times(&r2, m, y)));
        rat_clear(&a2); rat_clear(&r2);
    }
    Expr* e = surd_times(&c1, m, body);
    rat_clear(&a); rat_clear(&b); rat_clear(&A); rat_clear(&r); rat_clear(&c1); rat_clear(&z0); rat_clear(&mq);
    mpz_clear(m);
    return e;
}

// b·y split by the partial fractions of b·U: bU·y with den(bU) a power of U, the algebraic part
// G·y from Hermite steps at the other factors, and the simple poles r[i]/(d[i]·y), d[i] monic
typedef struct {
    QFn bU, G;
    Polynomial **d, **r;
    int count;
} PoleSplit;

static void pole_split_free(PoleSplit* S) {
    qfn_free(&S->bU); qfn_free(&S->G);
    for (int i = 0; i < S->count; i++) { free_polynomial(S->d[i]); free_polynomial(S->r[i]); }
    free(S->d); free(S->r);
}

// acc += p/q (takes p and q)
static void qfn_accumulate(QFn* acc, Polynomial* p, Polynomial* q) {
    QFn t = qfn_from(p, q), s = qfn_add(acc, &t, 1);
    qfn_free(acc); qfn_free(&t);
    *acc = s;
}

// a pole d^j of the integrand (d coprime to U)
typedef struct {
    const Polynomial *U, *d;
    int j;
} Pole;

// one Hermite step, j >= 2: p/(d^j·y) = D(Q·y/d^(j-1)) + p1/(d^(j-1)·y), where
// D(Q·y/d^(j-1)) = (Q'·U·d + Q·U'·d/2 - (j-1)·Q·d'·U)/(d^j·y) and Q ≡ -p/((j-1)·d'·U) mod d; returns p1
static Polynomial* hermite_step(const Pole* P, const Polynomial* p, Polynomial** Q) {
    Rational k, half;
    rat_init(&k); rat_init(&half);
    mpq_set_si(k.val, P->j - 1, 1);
    mpq_set_si(half.val, 1, 2);
    Polynomial *dd = poly_derivative(P->d), *dU = poly_derivative(P->U), *ddU = poly_multiply(dd, P->U);
    Polynomial *kddU = poly_scale(ddU, &k), *m = poly_mod(kddU, P->d);
    PolyBezout bz = poly_xgcd(m, P->d);                       // s·m ≡ 1 mod d
    Polynomial *ps = poly_multiply(p, bz.s), *q = poly_mod(ps, P->d), *z = create_polynomial();
    *Q = poly_subtract(z, q);
    Polynomial *dQ = poly_derivative(*Q), *a = poly_multiply(dQ, P->U), *ad = poly_multiply(a, P->d);
    Polynomial *b = poly_multiply(*Q, dU), *bd = poly_multiply(b, P->d), *bh = poly_scale(bd, &half);
    Polynomial *c = poly_multiply(*Q, kddU), *s1 = poly_add(ad, bh), *D = poly_subtract(s1, c);
    Polynomial *rest = poly_subtract(p, D), *p1 = poly_divide(rest, P->d);
    free_polynomial(dd); free_polynomial(dU); free_polynomial(ddU); free_polynomial(kddU); free_polynomial(m);
    free_polynomial(ps); free_polynomial(q); free_polynomial(z); free_polynomial(dQ); free_polynomial(a);
    free_polynomial(ad); free_polynomial(b); free_polynomial(bd); free_polynomial(bh); free_polynomial(c);
    free_polynomial(s1); free_polynomial(D); free_polynomial(rest);
    free_bezout(&bz);
    rat_clear(&k); rat_clear(&half);
    return p1;
}

// numerators[j-1]/F^j for a primitive factor F = l·d (d monic): numerators[j-1]/l^j over d^j
static Polynomial* over_monic(const PFTerm* T, int j, const Rational* l) {
    Rational lj;
    rat_init(&lj);
    rat_pow_int(&lj, l, -j);
    Polynomial* p = poly_scale(T->numerators[j - 1], &lj);
    rat_clear(&lj);
    return p;
}

// a factor F ≠ U of the denominator of b·U: Hermite steps from F^m down to F, then the simple pole
static void split_other(const SqCtx* c, const PFTerm* T, PoleSplit* S) {
    Rational l;
    rat_init(&l);
    poly_coeff(&l, T->factor, polynomial_degree(T->factor));
    Polynomial* d = poly_copy(T->factor);
    poly_make_monic(d);
    Polynomial* p = over_monic(T, T->power, &l);
    for (int j = T->power; j >= 2; j--) {
        Polynomial* Q;
        Polynomial *p1 = hermite_step(&(Pole){ c->U, d, j }, p, &Q), *nx = over_monic(T, j - 1, &l);
        qfn_accumulate(&S->G, Q, poly_pow(d, (unsigned)(j - 1)));
        free_polynomial(p);
        p = poly_add(p1, nx);
        free_polynomial(p1); free_polynomial(nx);
    }
    PolyDivResult qr = poly_divmod(p, d);                     // p/(d·y) = quo/y + rem/(d·y)
    qfn_accumulate(&S->bU, qr.quotient, poly_copy(c->U));
    S->d[S->count] = d;
    S->r[S->count++] = qr.remainder;
    free_polynomial(p);
    rat_clear(&l);
}

// the factor F of the partial fractions is U up to a constant: numerators[j-1]/F^j = n·κ^j/U^j, κ = U/F
static bool split_radicand(const SqCtx* c, const PFTerm* T, PoleSplit* S) {
    Rational lu, lf, k;
    rat_init(&lu); rat_init(&lf); rat_init(&k);
    poly_coeff(&lu, c->U, 2);
    poly_coeff(&lf, T->factor, polynomial_degree(T->factor));
    rat_div(&k, &lu, &lf);
    Polynomial* F = poly_scale(T->factor, &k);
    bool is_u = same_poly(F, c->U);
    for (int j = 1; is_u && j <= T->power; j++) {             // b gets n·κ^j/U^(j+1)
        Rational kj;
        rat_init(&kj);
        rat_pow_int(&kj, &k, j);
        qfn_accumulate(&S->bU, poly_scale(T->numerators[j - 1], &kj), poly_pow(c->U, (unsigned)(j + 1)));
        rat_clear(&kj);
    }
    free_polynomial(F);
    rat_clear(&lu); rat_clear(&lf); rat_clear(&k);
    return is_u;
}

static void split_poles(const SqCtx* c, const QFn* b, PoleSplit* S) {
    QFn u = qfn_poly(c->U), R = qfn_mul(b, &u);
    PartialFractionDecomp* pf = partial_fractions(R.num, R.den);
    *S = (PoleSplit){ qfn_zero(), qfn_zero(), malloc((size_t)(pf->num_terms + 1) * sizeof(Polynomial*)),
                      malloc((size_t)(pf->num_terms + 1) * sizeof(Polynomial*)), 0 };
    if (pf->poly_part) qfn_accumulate(&S->bU, poly_copy(pf->poly_part), poly_copy(c->U));
    for (int i = 0; i < pf->num_terms; i++)
        if (!split_radicand(c, &pf->terms[i], S)) split_other(c, &pf->terms[i], S);
    free_pf_decomp(pf);
    qfn_free(&u); qfn_free(&R);
}

// a reduction of b·y: the algebraic part alg = G·y, λ of λ·∫ 1/y, the simple poles, their log terms
typedef struct {
    PoleSplit S;
    SqFn alg;
    Rational lambda;
    Expr* logs;
} Reduction;

// the chain for b with den(b) a power of U: b·y = D(G·y) + λ/y with G = Σ Q_i/U^(m_i), into R
static bool reduce_powers(const SqCtx* c, const QFn* b, Reduction* R) {
    Polynomial* N = NULL;
    int m = 0;
    if (!as_reduction(c->U, b, &N, &m)) return false;
    QFn G0 = qfn_zero(), *G = &G0;
    Rational* lambda = &R->lambda;
    for (; m >= 1; m--) {
        Polynomial* Q;
        Polynomial* N1 = reduce_step(c->U, N, m, &Q);
        Polynomial* Um = poly_pow(c->U, (unsigned)m);
        QFn t = qfn_from(Q, Um), s = qfn_add(G, &t, 1);
        qfn_free(G); qfn_free(&t);
        *G = s;
        free_polynomial(N);
        N = N1;
    }
    Polynomial* Q0 = reduce_base(c->U, N, lambda);
    QFn q0 = qfn_poly(Q0), s = qfn_add(G, &q0, 1), all = qfn_add(&s, &R->S.G, 1);
    qfn_free(G); qfn_free(&q0); qfn_free(&s);
    R->alg = (SqFn){ qfn_zero(), all };                      // with the Hermite steps at the other poles
    free_polynomial(N); free_polynomial(Q0);
    return true;
}

// the exact check of the algebraic part: D(G·y) + λ/y + Σ r_i/(d_i·y) = b·y
static bool reduction_checks(const SqCtx* c, const Reduction* R, const QFn* b) {
    const PoleSplit* S = &R->S;
    SqFn d = sq_deriv(c, &R->alg);
    QFn acc = qfn_copy(&d.b);
    qfn_accumulate(&acc, poly_constant(&R->lambda), poly_copy(c->U));
    for (int i = 0; i < S->count; i++) qfn_accumulate(&acc, poly_copy(S->r[i]), poly_multiply(S->d[i], c->U));
    SqFn want = { qfn_zero(), qfn_copy(b) }, got = { qfn_copy(&d.a), acc };
    bool ok = sq_equal(&got, &want);
    sq_free(&d); sq_free(&want); sq_free(&got);
    return ok;
}

// F' ≈ f at points inside the real domain of y (central differences): a guard on the ln, atan and
// RootSum terms of the simple poles, which the exact check does not cover
static bool numeric_agrees(const Expr* F, const Expr* f, const Polynomial* U) {
    static const double u[] = { -0.93, -0.71, -0.52, -0.37, -0.18, 0.11, 0.29, 0.43, 0.61, 0.83 };
    Rational a, b, c0;
    rat_init(&a); rat_init(&b); rat_init(&c0);
    poly_coeff(&a, U, 2); poly_coeff(&b, U, 1); poly_coeff(&c0, U, 0);
    double A = rat_to_double(&a), B = rat_to_double(&b), C = rat_to_double(&c0);
    double mid = -B / (2 * A), disc = B * B - 4 * A * C, h = disc > 0 ? sqrt(disc) / (2 * fabs(A)) : 0;
    rat_clear(&a); rat_clear(&b); rat_clear(&c0);
    int checked = 0;
    for (size_t i = 0; i < 2 * (sizeof u / sizeof u[0]); i++) {
        double v = u[i % 10], x0 = A < 0 ? mid + h * v : mid + (h + 1) * (v + (v > 0 ? 1 : -1)) * (i < 10 ? 1 : 2);
        bool ok = true;
        double s = 1e-4 * (1 + fabs(x0)), g = expr_eval(f, x0, &ok);
        double d1 = (expr_eval(F, x0 + s, &ok) - expr_eval(F, x0 - s, &ok)) / (2 * s);
        double d2 = (expr_eval(F, x0 + s / 2, &ok) - expr_eval(F, x0 - s / 2, &ok)) / s;
        if (!ok || !isfinite(g) || fabs(g) > 1e6 || fabs(d1 - d2) > 1e-6 * (1 + fabs(g))) continue;
        if (fabs(d2 - g) > 1e-6 * (1 + fabs(g))) return false;
        checked++;
    }
    return checked >= 2;
}

// the ln/atan/RootSum terms of the simple poles r_i/(d_i·y) in R->logs (NULL when there are none);
// false with *why when one is unsupported
static bool simple_poles(const SqCtx* c, Reduction* R, const char** why) {
    const PoleSplit* S = &R->S;
    Expr** out = &R->logs;
    *out = NULL;
    for (int i = 0; i < S->count; i++) {
        if (poly_is_zero(S->r[i])) continue;
        Expr* t = sqrt_log_part(&(SqrtLog){ c->U, S->r[i], S->d[i] }, why);
        if (!t) { free_expr(*out); *out = NULL; return false; }
        *out = expr_append_sum(*out, t);
    }
    return true;
}

// ∫ b·y on a conic without a rational point: Hermite steps at every pole (U's powers by the
// chain to ∫ 1/y), λ·∫ 1/y, and the simple poles' logarithmic part; the algebraic part is checked
// exactly, the whole answer numerically; NULL with *why otherwise
static Expr* integrate_by_reduction(const SqCtx* c, const QFn* b, const char** why) {
    Reduction R = { .logs = NULL };
    split_poles(c, b, &R.S);
    rat_init(&R.lambda);
    Expr* F = NULL;
    if (!reduce_powers(c, &R.S.bU, &R)) {
        *why = "internal error: the square-root reduction left a denominator that is not a power of U";
        pole_split_free(&R.S); rat_clear(&R.lambda);
        return NULL;
    }
    bool ok = reduction_checks(c, &R, b);
    if (!ok) *why = "internal check failed: the square-root reduction does not differentiate back";
    if (ok && simple_poles(c, &R, why)) {
        F = sq_is_zero(&R.alg) ? NULL : sq_expr(&R.alg, c->U);
        if (!rat_is_zero(&R.lambda)) F = expr_append_sum(F, base_integral(c->U, &R.lambda));
        if (R.logs) F = expr_append_sum(F, R.logs);
        if (!F) F = create_const(0);
        SqFn fb = { qfn_zero(), qfn_copy(b) };
        Expr* f = sq_expr(&fb, c->U);
        if (R.S.count && !numeric_agrees(F, f, c->U)) {
            *why = "internal check failed: the logarithmic part of the square-root integral does not differentiate back";
            free_expr(F);
            F = NULL;
        }
        free_expr(f);
        sq_free(&fb);
    }
    sq_free(&R.alg);
    pole_split_free(&R.S);
    rat_clear(&R.lambda);
    return F;
}

// ∫ b·y for deg U = 1 or 2: R(t) = b(X)·Y·X', integrated over Q, then t = T(x, y)
static Expr* integrate_root_part(const SqCtx* c, const QFn* b, const char** why) {
    Param P;
    if (!param_of(c->U, &P)) return integrate_by_reduction(c, b, why);    // no rational point: reduce
    Expr* F = NULL;
    if (!param_inverts(c, &P)) *why = "internal check failed: the parametrization of the square root does not invert";
    else {
        QFn bX = qfn_at_q(b, &P.X), dX = qfn_derivative(&P.X), bY = qfn_mul(&bX, &P.Y), R = qfn_mul(&bY, &dX);
        Expr* r = qfn_expr(&R);
        Expr* Ft = integrate(r, "x");
        if (!Ft) *why = integrate_error();
        else { F = back_sum(Ft, &(Back){ c, &P.T }); free_expr(Ft); }
        free_expr(r);
        qfn_free(&bX); qfn_free(&dX); qfn_free(&bY); qfn_free(&R);
    }
    param_free(&P);
    return F;
}

// ∫ a for a ∈ Q(x), NULL with *why on failure
static Expr* integrate_rational_part(const QFn* a, const char** why) {
    Expr* e = qfn_expr(a);
    Expr* F = integrate(e, "x");
    if (!F) *why = integrate_error();
    free_expr(e);
    return F;
}

// ∫ b·sqrt(m) = sqrt(m)·∫ b for a constant radicand (a rational ∫ b written as q·sqrt(m))
static Expr* integrate_constant_root(const SqCtx* c, const QFn* b, const char** why) {
    Expr* F = integrate_rational_part(b, why);
    if (!F) return NULL;
    Expr* k = canonicalize(F);
    Expr* r = NULL;
    if (k && (k->type == EXPR_POLY || k->type == EXPR_RATIONAL_FN)) {
        QFn q = k->type == EXPR_POLY ? qfn_poly(k->poly) : qfn_from(poly_copy(k->ratfn.num), poly_copy(k->ratfn.den));
        r = times_root(&q, c->U);
        qfn_free(&q);
    } else r = create_mul(create_sqrt(create_poly_expr(poly_copy(c->U))), expr_copy(F));
    free_expr(k); free_expr(F);
    return r;
}

RischResult integrate_sqrt(const Expr* f) {
    SqCtx c;
    SqFn g;
    RischResult r = { RISCH_UNSUPPORTED, NULL, NULL };
    if (!form_of(f, &c, &g)) { r.reason = c.why; free_polynomial(c.U); return r; }
    const char* why = NULL;
    int deg = c.U ? polynomial_degree(c.U) : -1;
    Expr *A = qfn_is_zero(&g.a) ? NULL : integrate_rational_part(&g.a, &why), *Bp = NULL;
    bool ok = qfn_is_zero(&g.a) || A;
    if (ok && !qfn_is_zero(&g.b)) {
        if (deg > 2) why = "square roots of polynomials of degree 3 or more (elliptic integrals and beyond) are not supported";
        else Bp = deg == 0 ? integrate_constant_root(&c, &g.b, &why) : integrate_root_part(&c, &g.b, &why);
        ok = Bp != NULL;
    }
    if (ok) {
        r.status = RISCH_OK;
        r.result = A && Bp ? create_add(A, Bp) : A ? A : Bp ? Bp : create_const(0);
    } else { free_expr(A); free_expr(Bp); r.reason = why; }
    sq_free(&g);
    free_polynomial(c.U);
    return r;
}

// --- values at a point, for definite integrals ---

bool sqrt_shape(const Expr* f, SqrtShape* s, const char** why) {
    SqCtx c;
    SqFn g;
    *s = (SqrtShape){ NULL, NULL, NULL, NULL, NULL };
    if (!form_of(f, &c, &g)) { *why = c.why; free_polynomial(c.U); return false; }
    *s = (SqrtShape){ c.U, poly_copy(g.a.num), poly_copy(g.a.den), poly_copy(g.b.num), poly_copy(g.b.den) };
    sq_free(&g);
    return true;
}

void sqrt_shape_free(SqrtShape* s) {
    free_polynomial(s->U); free_polynomial(s->num_a); free_polynomial(s->den_a);
    free_polynomial(s->num_b); free_polynomial(s->den_b);
    *s = (SqrtShape){ NULL, NULL, NULL, NULL, NULL };
}

// p = (x - x0)^k·q with q(x0) ≠ 0: returns q, sets *k (p ≠ 0)
static Polynomial* strip_root(const Polynomial* p, const Rational* x0, int* k) {
    Rational one = rat_one_val(), m;
    rat_init(&m);
    rat_neg(&m, x0);
    Polynomial *lin = create_polynomial(), *q = poly_copy(p);
    add_term(lin, &one, 1); add_term(lin, &m, 0);
    for (*k = 0;; (*k)++) {
        PolyDivResult d = poly_divmod(q, lin);
        bool divides = poly_is_zero(d.remainder);
        free_polynomial(d.remainder);
        if (!divides) { free_polynomial(d.quotient); break; }
        free_polynomial(q);
        q = d.quotient;
    }
    free_polynomial(lin);
    rat_clear(&one); rat_clear(&m);
    return q;
}

// y0 = sqrt(u0) for u0 >= 0 in Q(sqrt(m)): the context (U0 = m, or none) and y0 as an element
static SqFn root_constant(const Rational* u0, SqCtx* c0) {
    *c0 = (SqCtx){ NULL, NULL };
    if (rat_is_zero(u0)) return sq_zero();
    Polynomial* p = poly_constant(u0);
    SqFn y0;
    radical(c0, p, &y0);                                      // rational, or r·sqrt(m) with c0->U = m
    free_polynomial(p);
    return y0;
}

// the element A0 + B0·y0 of Q(sqrt(u0)) for polynomials A, B at x0
static SqFn at_point(const SqCtx* c0, const SqFn* y0, const Polynomial* A, const Polynomial* B, const Rational* x0) {
    Rational va, vb;
    rat_init(&va); rat_init(&vb);
    poly_eval_rat(&va, A, x0); poly_eval_rat(&vb, B, x0);
    SqFn a = sq_rat(qfn_rat(&va)), b = sq_rat(qfn_rat(&vb)), by = sq_mul(c0, &b, y0), r = sq_add(&a, &by, 1);
    sq_free(&a); sq_free(&b); sq_free(&by);
    rat_clear(&va); rat_clear(&vb);
    return r;
}

// g = (A + B·y)/D at x0: directly when D(x0) ≠ 0, else (A² - B²·U)/(D·(A - B·y)) with the factor
// (x - x0)^k of D cancelled from A² - B²·U; false when g is infinite at x0
static bool field_value(const SqCtx* c, const SqFn* g, const Rational* x0, Expr** out) {
    Polynomial* D = poly_lcm(g->a.den, g->b.den);
    QFn d = qfn_poly(D), Aq = qfn_mul(&g->a, &d), Bq = qfn_mul(&g->b, &d);
    Rational u0, dv;
    rat_init(&u0); rat_init(&dv);
    if (c->U) poly_eval_rat(&u0, c->U, x0);
    SqCtx c0;
    SqFn y0 = root_constant(&u0, &c0);
    poly_eval_rat(&dv, D, x0);
    bool ok = true;
    SqFn v;
    if (!rat_is_zero(&dv)) {                                  // (A0 + B0·y0)/D0
        SqFn n = at_point(&c0, &y0, Aq.num, Bq.num, x0), id = sq_rat(qfn_rat(&dv)), inv = sq_inv(&c0, &id);
        v = sq_mul(&c0, &n, &inv);
        sq_free(&n); sq_free(&id); sq_free(&inv);
    } else {
        Polynomial *A2 = poly_multiply(Aq.num, Aq.num), *B2 = poly_multiply(Bq.num, Bq.num);
        Polynomial *B2U = c->U ? poly_multiply(B2, c->U) : create_polynomial(), *N = poly_subtract(A2, B2U);
        int kd = 0, kn = 0;
        Polynomial *D1 = strip_root(D, x0, &kd), *N1 = poly_is_zero(N) ? poly_copy(N) : strip_root(N, x0, &kn);
        ok = !poly_is_zero(N) && kn >= kd;
        if (ok) {                                             // N1·(x - x0)^(kn - kd) / (D1·(A - B·y)) at x0
            Rational nv, d1;
            rat_init(&nv); rat_init(&d1);
            poly_eval_rat(&nv, N1, x0);
            if (kn > kd) mpq_set_ui(nv.val, 0, 1);
            poly_eval_rat(&d1, D1, x0);
            Rational m1 = rat_from_int(-1);
            Polynomial* nB = poly_scale(Bq.num, &m1);
            SqFn conj = at_point(&c0, &y0, Aq.num, nB, x0), dd = sq_rat(qfn_rat(&d1)), den = sq_mul(&c0, &dd, &conj);
            ok = !sq_is_zero(&den);
            if (ok) {
                SqFn inv = sq_inv(&c0, &den), nn = sq_rat(qfn_rat(&nv));
                v = sq_mul(&c0, &nn, &inv);
                sq_free(&inv); sq_free(&nn);
            }
            sq_free(&conj); sq_free(&dd); sq_free(&den);
            free_polynomial(nB);
            rat_clear(&nv); rat_clear(&d1); rat_clear(&m1);
        }
        free_polynomial(A2); free_polynomial(B2); free_polynomial(B2U); free_polynomial(N); free_polynomial(D1); free_polynomial(N1);
    }
    if (ok) { *out = c0.U ? sq_expr(&v, c0.U) : qfn_expr(&v.a); sq_free(&v); }
    sq_free(&y0);
    free_polynomial(c0.U); free_polynomial(D);
    qfn_free(&d); qfn_free(&Aq); qfn_free(&Bq);
    rat_clear(&u0); rat_clear(&dv);
    return ok;
}

Expr* sqrt_value_at(const Expr* g, const Polynomial* U, const Rational* x0) {
    SqCtx c = { U ? poly_copy(U) : NULL, NULL };
    SqFn f;
    Expr* out = NULL;
    if (conv(g, &c, &f)) {
        if (!field_value(&c, &f, x0, &out)) out = NULL;
        sq_free(&f);
    }
    free_polynomial(c.U);
    return out;
}
