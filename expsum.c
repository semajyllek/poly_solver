// Exponential sums over Q(i)(x): conversion from expressions, integration term by term with
// the Risch differential equation over Q(i), and the real form for output.
//
// Two exponents that differ by a nonzero constant c (e^{i(x+1)} and e^{ix}) stay separate terms.
// That keeps proofs exact: e^c is transcendental (Lindemann), and the image of y ↦ y' + w'·y is
// defined over Q(i), so c₁ + e^c·c₂ lies in it only if c₁ and c₂ both do.
#include "expsum.h"
#include "cfield.h"
#include <stdlib.h>

#define MAX_POWER 256

typedef struct {
    CPoly w;
    CFn c;
} ETerm;

typedef struct {
    ETerm* t;
    int n;
} ESum;

static ESum es_empty(void) { return (ESum){ NULL, 0 }; }

static void es_free(ESum* s) {
    for (int i = 0; i < s->n; i++) { cp_free(&s->t[i].w); cfn_free(&s->t[i].c); }
    free(s->t);
    *s = es_empty();
}

// s += c·e^w; takes w and c
static void es_add_term(ESum* s, CPoly w, CFn c) {
    for (int i = 0; i < s->n; i++) {
        if (!cp_equal(&s->t[i].w, &w)) continue;
        CFn sum = cfn_add(&s->t[i].c, &c, 1);
        cfn_free(&s->t[i].c); cfn_free(&c); cp_free(&w);
        s->t[i].c = sum;
        if (cfn_is_zero(&sum)) {                          // drop it: the last term takes its place
            cp_free(&s->t[i].w); cfn_free(&s->t[i].c);
            s->t[i] = s->t[--s->n];
        }
        return;
    }
    if (cfn_is_zero(&c)) { cp_free(&w); cfn_free(&c); return; }
    s->t = realloc(s->t, (size_t)(s->n + 1) * sizeof(ETerm));
    s->t[s->n++] = (ETerm){ w, c };
}

static ESum es_term(CPoly w, CFn c) {
    ESum s = es_empty();
    es_add_term(&s, w, c);
    return s;
}

static ESum es_add(const ESum* a, const ESum* b, int sign) {
    ESum r = es_empty();
    for (int i = 0; i < a->n; i++) es_add_term(&r, cp_copy(&a->t[i].w), cfn_copy(&a->t[i].c));
    CFn m = cfn_const(sign, 0, 1);
    for (int i = 0; i < b->n; i++) es_add_term(&r, cp_copy(&b->t[i].w), cfn_mul(&b->t[i].c, &m));
    cfn_free(&m);
    return r;
}

// products add exponents
static ESum es_mul(const ESum* a, const ESum* b) {
    ESum r = es_empty();
    for (int i = 0; i < a->n; i++)
        for (int j = 0; j < b->n; j++) es_add_term(&r, cp_add(&a->t[i].w, &b->t[j].w, 1), cfn_mul(&a->t[i].c, &b->t[j].c));
    return r;
}

// D(c·e^w) = (c' + w'·c)·e^w
static ESum es_derivative(const ESum* s) {
    ESum r = es_empty();
    for (int i = 0; i < s->n; i++) {
        CPoly dw = cp_derivative(&s->t[i].w);
        CFn wq = cfn_poly(&dw), wc = cfn_mul(&wq, &s->t[i].c), dc = cfn_derivative(&s->t[i].c);
        es_add_term(&r, cp_copy(&s->t[i].w), cfn_add(&dc, &wc, 1));
        cp_free(&dw); cfn_free(&wq); cfn_free(&wc); cfn_free(&dc);
    }
    return r;
}

static bool es_equal(const ESum* a, const ESum* b) {
    ESum d = es_add(a, b, -1);
    bool z = d.n == 0;
    es_free(&d);
    return z;
}

// --- expression -> exponential sum ---

typedef struct {
    const char* unsupported;
} EConv;

static ESum conv(const Expr* e, EConv* C);

// the real polynomial an exp/sin/cos argument must be
static Polynomial* poly_arg(const Expr* arg, EConv* C) {
    Expr* c = canonicalize((Expr*)arg);                   // canonicalize does not modify arg
    Polynomial* p = c && c->type == EXPR_POLY ? poly_copy(c->poly) : NULL;
    free_expr(c);
    if (!p) C->unsupported = "exp, sin or cos of a non-polynomial argument (such as sin(1/x) or exp(exp(x))) is not supported yet";
    return p;
}

// exp(u) = e^u; cos(u) = (e^{iu} + e^{-iu})/2; sin(u) = -i/2·e^{iu} + i/2·e^{-iu}
static ESum conv_fn(const Expr* e, EConv* C) {
    Polynomial* u = poly_arg(e->fn_arg, C);
    if (!u) return es_empty();
    ESum s = es_empty();
    if (e->type == EXPR_EXP) es_add_term(&s, (CPoly){ u, create_polynomial() }, cfn_const(1, 0, 1));
    else {
        CPoly w = { create_polynomial(), u }, mw = cp_neg(&w);
        bool sine = e->type == EXPR_SIN;
        es_add_term(&s, w, sine ? cfn_const(0, -1, 2) : cfn_const(1, 0, 2));
        es_add_term(&s, mw, sine ? cfn_const(0, 1, 2) : cfn_const(1, 0, 2));
    }
    return s;
}

// 1/(c·e^w) = (1/c)·e^{-w}: only single terms can be divisors
static ESum es_inverse(const ESum* s, EConv* C) {
    if (s->n == 0) { C->unsupported = "division by zero"; return es_empty(); }
    if (s->n > 1) {
        C->unsupported = "division by a sum involving sin, cos or several exponentials (such as 1/(1 + cos(x)) or tan(x)) is not supported yet";
        return es_empty();
    }
    return es_term(cp_neg(&s->t[0].w), cfn_inv(&s->t[0].c));
}

// a small integer exponent, or false
static bool small_integer(const Expr* e, long* k) {
    Expr* n = canonicalize((Expr*)e);
    bool ok = n && n->type == EXPR_POLY && polynomial_degree(n->poly) <= 0;
    if (ok) {
        Rational c;
        rat_init(&c);
        poly_coeff(&c, n->poly, 0);
        ok = mpz_cmp_ui(mpq_denref(c.val), 1) == 0 && mpz_cmpabs_ui(mpq_numref(c.val), MAX_POWER) <= 0;
        if (ok) *k = mpz_get_si(mpq_numref(c.val));
        rat_clear(&c);
    }
    free_expr(n);
    return ok;
}

static ESum conv_pow(const Expr* e, EConv* C) {
    long k = 0;
    if (!small_integer(e->binop.right, &k)) { C->unsupported = "exponents must be small integers"; return es_empty(); }
    ESum b = conv(e->binop.left, C);
    if (k < 0) { ESum inv = es_inverse(&b, C); es_free(&b); b = inv; k = -k; }
    ESum r = es_term(cp_zero(), cfn_const(1, 0, 1));
    for (long i = 0; i < k && !C->unsupported; i++) { ESum p = es_mul(&r, &b); es_free(&r); r = p; }
    es_free(&b);
    return r;
}

static ESum conv_binop(const Expr* e, EConv* C) {
    if (e->binop.op == OP_POW) return conv_pow(e, C);
    ESum l = conv(e->binop.left, C), r = conv(e->binop.right, C), out = es_empty();
    if (!C->unsupported) switch (e->binop.op) {
        case OP_ADD: out = es_add(&l, &r, 1); break;
        case OP_SUB: out = es_add(&l, &r, -1); break;
        case OP_MUL: out = es_mul(&l, &r); break;
        default: { ESum inv = es_inverse(&r, C); out = es_mul(&l, &inv); es_free(&inv); }
    }
    es_free(&l); es_free(&r);
    return out;
}

static ESum conv(const Expr* e, EConv* C) {
    if (C->unsupported) return es_empty();
    if (expr_has_ln(e)) { C->unsupported = "ln together with sin/cos or several exponentials is not supported yet"; return es_empty(); }
    if (!expr_has_exp(e) && !expr_has_trig(e)) {          // a rational function of x
        Expr* c = canonicalize((Expr*)e);
        if (!c) { C->unsupported = "not a rational function of x, exp, sin and cos"; return es_empty(); }
        CPoly n = { poly_copy(c->type == EXPR_POLY ? c->poly : c->ratfn.num), create_polynomial() };
        CFn a = c->type == EXPR_POLY ? cfn_poly(&n) : cfn_real(c->ratfn.num, c->ratfn.den);
        cp_free(&n);
        free_expr(c);
        return es_term(cp_zero(), a);
    }
    switch (e->type) {
        case EXPR_BINOP: return conv_binop(e, C);
        case EXPR_NEG: { ESum z = es_empty(), f = conv(e->operand, C), r = es_add(&z, &f, -1); es_free(&f); return r; }
        case EXPR_EXP: case EXPR_SIN: case EXPR_COS: return conv_fn(e, C);
        default: C->unsupported = "only +, -, *, /, integer powers, exp, sin and cos are supported here"; return es_empty();
    }
}

// a nonzero constant exponent is a transcendental constant factor (exp(3), cos(1))
static bool constant_exponent(const ESum* s) {
    for (int i = 0; i < s->n; i++)
        if (cp_deg(&s->t[i].w) == 0) return true;
    return false;
}

static ESum es_from(const Expr* e, const char** why) {
    EConv C = { NULL };
    ESum s = conv(e, &C);
    if (!C.unsupported && constant_exponent(&s))
        C.unsupported = "exp, sin or cos of a constant (a transcendental number such as exp(3) or cos(1)) is not supported yet";
    if (C.unsupported) es_free(&s);
    *why = C.unsupported;
    return s;
}

// --- the real form: c·e^{E+iF} + d·e^{E-iF} = e^E·((c + d)·cos F + i(c - d)·sin F) ---

// the sign of p's leading coefficient (0 for p = 0)
static int lead_sign(const Polynomial* p) {
    if (poly_is_zero(p)) return 0;
    return rat_is_negative(&p->head->coeff) ? -1 : 1;
}

// a canonical order on polynomials: degree, then coefficients from the top
static int poly_order(const Polynomial* a, const Polynomial* b) {
    int da = poly_is_zero(a) ? -1 : polynomial_degree(a), db = poly_is_zero(b) ? -1 : polynomial_degree(b);
    if (da != db) return da < db ? -1 : 1;
    for (const Term *s = a->head, *t = b->head; s && t; s = s->next, t = t->next) {
        if (s->exponent != t->exponent) return s->exponent < t->exponent ? -1 : 1;
        int c = rat_cmp(&s->coeff, &t->coeff);
        if (c) return c;
    }
    return 0;
}

// the terms with exponents E ± i·F, F with a positive leading coefficient (or 0)
typedef struct {
    Polynomial *E, *F;
    CFn c, d;     // coefficients of e^{E+iF} and e^{E-iF}
} Pair;

static int pair_order(const void* x, const void* y) {
    const Pair *a = x, *b = y;
    int f = poly_order(a->F, b->F);
    return f ? f : poly_order(a->E, b->E);
}

typedef struct {
    Pair* p;
    int n;
} Pairs;

static Pairs pairs_of(const ESum* s) {
    Pairs P = { malloc((size_t)(s->n > 0 ? s->n : 1) * sizeof(Pair)), 0 };
    for (int i = 0; i < s->n; i++) {
        int sg = lead_sign(s->t[i].w.im);
        CPoly key = sg < 0 ? cp_neg(&s->t[i].w) : cp_copy(&s->t[i].w);
        if (sg < 0) { Polynomial* E = poly_copy(s->t[i].w.re); free_polynomial(key.re); key.re = E; } // conjugate: E stays
        int j = 0;
        while (j < P.n && !(poly_order(P.p[j].E, key.re) == 0 && poly_order(P.p[j].F, key.im) == 0)) j++;
        if (j == P.n) P.p[P.n++] = (Pair){ poly_copy(key.re), poly_copy(key.im), cfn_zero(), cfn_zero() };
        CFn* slot = sg < 0 ? &P.p[j].d : &P.p[j].c;
        cfn_free(slot);
        *slot = cfn_copy(&s->t[i].c);
        cp_free(&key);
    }
    qsort(P.p, (size_t)P.n, sizeof(Pair), pair_order);
    return P;
}

static void pairs_free(Pairs* P) {
    for (int i = 0; i < P->n; i++) { free_polynomial(P->p[i].E); free_polynomial(P->p[i].F); cfn_free(&P->p[i].c); cfn_free(&P->p[i].d); }
    free(P->p);
}

// a constant polynomial equal to ±1
static int unit_sign(const Polynomial* p) {
    if (poly_is_zero(p) || polynomial_degree(p) != 0) return 0;
    const Rational* c = &p->head->coeff;
    if (mpz_cmpabs_ui(mpq_numref(c->val), 1) || mpz_cmp_ui(mpq_denref(c->val), 1)) return 0;
    return rat_is_negative(c) ? -1 : 1;
}

// (num/den)·core, sign outside for a coefficient of -1; takes core; NULL when num = 0
static Expr* scaled(const Polynomial* num, const Polynomial* den, Expr* core) {
    if (poly_is_zero(num)) { free_expr(core); return NULL; }
    int u = unit_sign(num);
    Expr* e = u ? core : create_mul(create_poly_expr(poly_copy(num)), core);
    if (polynomial_degree(den) > 0) e = create_div(e, create_poly_expr(poly_copy(den)));
    return u < 0 ? create_neg(e) : e;
}

// exp(E)·trig(F), each factor present only when nonzero; trig NULL for none
static Expr* core_expr(const Pair* p, Expr* (*trig)(Expr*)) {
    Expr* ex = poly_is_zero(p->E) ? NULL : create_exp(create_poly_expr(poly_copy(p->E)));
    Expr* tr = trig ? trig(create_poly_expr(poly_copy(p->F))) : NULL;
    if (!ex) return tr ? tr : create_const(1);
    return tr ? create_mul(ex, tr) : ex;
}

static Expr* plus(Expr* acc, Expr* t) { return t ? expr_append_sum(acc, t) : acc; }

// e^E·(sin coefficient·sin F + cos coefficient·cos F); false when the pair is not real
static bool pair_expr(const Pair* p, Expr** acc) {
    if (poly_is_zero(p->F)) {                             // c·e^E, c real
        if (!poly_is_zero(p->c.n.im)) return false;
        *acc = plus(*acc, scaled(p->c.n.re, p->c.den, core_expr(p, NULL)));
        return true;
    }
    CFn sum = cfn_add(&p->c, &p->d, 1), dif = cfn_add(&p->c, &p->d, -1);
    bool real = poly_is_zero(sum.n.im) && poly_is_zero(dif.n.re);
    if (real) {                                           // i·(c - d) = -Im(c - d) for real pairs
        Rational m1 = rat_from_int(-1);
        Polynomial* sn = poly_scale(dif.n.im, &m1);
        *acc = plus(*acc, scaled(sn, dif.den, core_expr(p, create_sin)));
        *acc = plus(*acc, scaled(sum.n.re, sum.den, core_expr(p, create_cos)));
        free_polynomial(sn);
        rat_clear(&m1);
    }
    cfn_free(&sum); cfn_free(&dif);
    return real;
}

// s as a real expression (0 when empty); NULL when s is not a real function
static Expr* es_real_expr(const ESum* s) {
    Pairs P = pairs_of(s);
    Expr* acc = NULL;
    bool ok = true;
    for (int i = 0; ok && i < P.n; i++) ok = pair_expr(&P.p[i], &acc);
    pairs_free(&P);
    if (!ok) { free_expr(acc); return NULL; }
    return acc ? acc : create_const(0);
}

// --- integration ---

static RischResult failure(RischStatus s, const char* why) { return (RischResult){ s, NULL, why }; }

// the terms of s with exponent w or its conjugate (the real part a proof is about)
static ESum conjugate_pair(const ESum* s, const CPoly* w) {
    CPoly cw = { poly_copy(w->re), NULL };
    Rational m1 = rat_from_int(-1);
    cw.im = poly_scale(w->im, &m1);
    ESum r = es_empty();
    for (int i = 0; i < s->n; i++)
        if (cp_equal(&s->t[i].w, w) || cp_equal(&s->t[i].w, &cw)) es_add_term(&r, cp_copy(&s->t[i].w), cfn_copy(&s->t[i].c));
    cp_free(&cw);
    rat_clear(&m1);
    return r;
}

// Y = Σ y_w·e^w over w ≠ 0 with y_w' + w'·y_w = c_w; false (with *bad the failing index) if one has none
static bool solve_terms(const ESum* s, ESum* Y, int* bad) {
    for (int i = 0; i < s->n; i++) {
        if (cp_is_zero(&s->t[i].w)) continue;
        CPoly dw = cp_derivative(&s->t[i].w);
        CFn y;
        bool ok = c_solve_rde(&dw, &s->t[i].c, &y);
        cp_free(&dw);
        if (!ok) { *bad = i; return false; }
        es_add_term(Y, cp_copy(&s->t[i].w), y);
    }
    return true;
}

// the e^0 coefficient (a real rational function) as an expression
static Expr* rational_part(const ESum* s) {
    for (int i = 0; i < s->n; i++) {
        if (!cp_is_zero(&s->t[i].w)) continue;
        if (polynomial_degree(s->t[i].c.den) == 0) return create_poly_expr(poly_copy(s->t[i].c.n.re));
        return create_div(create_poly_expr(poly_copy(s->t[i].c.n.re)), create_poly_expr(poly_copy(s->t[i].c.den)));
    }
    return create_const(0);
}

static bool is_zero_expr(const Expr* e) {
    return (e->type == EXPR_CONST && rat_is_zero(&e->constant)) || (e->type == EXPR_POLY && poly_is_zero(e->poly));
}

// ∫ a_0 + Σ y_w·e^w in real form, after an exact check D(Σ y_w·e^w) = f - a_0
static RischResult assemble(const ESum* s, const ESum* Y) {
    ESum dY = es_derivative(Y), rest = es_empty(), a0 = es_empty();
    for (int i = 0; i < s->n; i++) es_add_term(cp_is_zero(&s->t[i].w) ? &a0 : &rest, cp_copy(&s->t[i].w), cfn_copy(&s->t[i].c));
    bool same = es_equal(&dY, &rest), real = a0.n == 0 || poly_is_zero(a0.t[0].c.n.im);
    es_free(&dY); es_free(&rest); es_free(&a0);
    if (!same || !real) return failure(RISCH_UNSUPPORTED, "internal check failed: the exp/sin/cos part does not differentiate back");
    Expr* r0 = rational_part(s);
    Expr* I0 = integrate(r0, "x");                       // self-checked
    free_expr(r0);
    if (!I0) return failure(RISCH_UNSUPPORTED, integrate_error());
    Expr* ey = es_real_expr(Y);
    if (!ey) { free_expr(I0); return failure(RISCH_UNSUPPORTED, "internal check failed: the antiderivative is not real"); }
    Expr* r = is_zero_expr(ey) ? (free_expr(ey), I0) : is_zero_expr(I0) ? (free_expr(I0), ey) : expr_append_sum(I0, ey);
    return (RischResult){ RISCH_OK, r, NULL };
}

RischResult integrate_expsum(const Expr* f) {
    const char* why = NULL;
    ESum s = es_from(f, &why);
    if (why) return failure(RISCH_UNSUPPORTED, why);
    ESum Y = es_empty();
    int bad = 0;
    RischResult r;
    if (solve_terms(&s, &Y, &bad)) r = assemble(&s, &Y);
    else {
        ESum part = conjugate_pair(&s, &s.t[bad].w);
        r = (RischResult){ RISCH_NONELEMENTARY, es_real_expr(&part),
                           "has no elementary antiderivative (Risch: y' + w'·y = c has no rational solution over Q(i))" };
        es_free(&part);
        if (!r.result) r = failure(RISCH_UNSUPPORTED, "the integrand is not a real function");
    }
    es_free(&s); es_free(&Y);
    return r;
}

Expr* differentiate_expsum(const Expr* f, const char** reason) {
    ESum s = es_from(f, reason);
    if (*reason) return NULL;
    ESum d = es_derivative(&s);
    Expr* r = es_real_expr(&d);
    if (!r) *reason = "the expression is not a real function";
    es_free(&s); es_free(&d);
    return r;
}

Expr* simplify_expsum(const Expr* f, const char** reason) {
    ESum s = es_from(f, reason);
    if (*reason) return NULL;
    Expr* r = es_real_expr(&s);
    if (!r) *reason = "the expression is not a real function";
    es_free(&s);
    return r;
}
