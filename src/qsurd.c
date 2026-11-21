// Arithmetic in Q(sqrt(m)), m a positive integer that is not a square.
#include "qsurd.h"
#include <stdlib.h>

void qs_init(QSurd* a) { rat_init(&a->p); rat_init(&a->q); a->m = 0; }
void qs_clear(QSurd* a) { rat_clear(&a->p); rat_clear(&a->q); }

void qs_set(QSurd* r, const QSurd* a) { rat_set(&r->p, &a->p); rat_set(&r->q, &a->q); r->m = a->m; }

// m = 0 whenever q = 0, so a rational fits every field
static void normalize(QSurd* a) { if (rat_is_zero(&a->q)) a->m = 0; }

// the common radicand of a and b, or -1 when they lie in different fields
static long field_of(const QSurd* a, const QSurd* b) {
    if (!a->m || !b->m || a->m == b->m) return a->m ? a->m : b->m;
    return -1;
}

// --- parsing ---

// a single term c·sqrt(m)^k (k = 0 or 1), from constants, one sqrt, products and quotients
static bool parse_term(const Expr* e, Rational* c, long* m) {
    if (e->type == EXPR_CONST) { rat_set(c, &e->constant); *m = 0; return true; }
    if (e->type == EXPR_NEG) { bool ok = parse_term(e->operand, c, m); rat_neg(c, c); return ok; }
    if (e->type == EXPR_SQRT) {
        const Expr* r = e->fn_arg;
        bool integer = r->type == EXPR_CONST && mpz_cmp_ui(mpq_denref(r->constant.val), 1) == 0
                    && mpz_sgn(mpq_numref(r->constant.val)) > 0 && mpz_fits_slong_p(mpq_numref(r->constant.val));
        if (!integer) return false;
        rat_set_si(c, 1, 1);
        *m = mpz_get_si(mpq_numref(r->constant.val));
        return true;
    }
    if (e->type != EXPR_BINOP || (e->binop.op != OP_MUL && e->binop.op != OP_DIV)) return false;
    Rational cr;
    rat_init(&cr);
    long mr;
    bool ok = parse_term(e->binop.left, c, m) && parse_term(e->binop.right, &cr, &mr);
    if (ok && e->binop.op == OP_MUL) ok = !(*m && mr);   // at most one sqrt factor
    if (ok && e->binop.op == OP_DIV) ok = !mr && !rat_is_zero(&cr);
    if (ok) {
        if (e->binop.op == OP_MUL) { rat_mul(c, c, &cr); if (mr) *m = mr; }
        else rat_div(c, c, &cr);
    }
    rat_clear(&cr);
    return ok;
}

// out += sign·e
static bool parse_sum(const Expr* e, int sign, QSurd* out) {
    bool add = e->type == EXPR_BINOP && e->binop.op == OP_ADD, sub = e->type == EXPR_BINOP && e->binop.op == OP_SUB;
    if (add || sub) return parse_sum(e->binop.left, sign, out) && parse_sum(e->binop.right, sub ? -sign : sign, out);
    Rational c;
    rat_init(&c);
    long m;
    bool ok = parse_term(e, &c, &m) && (!m || !out->m || m == out->m);
    if (ok) {
        if (sign < 0) rat_neg(&c, &c);
        if (m) { rat_add(&out->q, &out->q, &c); out->m = m; }
        else rat_add(&out->p, &out->p, &c);
    }
    rat_clear(&c);
    return ok;
}

bool qs_parse(const Expr* e, QSurd* out) {
    rat_set_si(&out->p, 0, 1); rat_set_si(&out->q, 0, 1); out->m = 0;
    bool ok = parse_sum(e, 1, out);
    normalize(out);
    return ok;
}

// --- arithmetic ---

static bool add_signed(QSurd* r, const QSurd* a, const QSurd* b, int sign) {
    long m = field_of(a, b);
    if (m < 0) return false;
    if (sign > 0) { rat_add(&r->p, &a->p, &b->p); rat_add(&r->q, &a->q, &b->q); }
    else { rat_sub(&r->p, &a->p, &b->p); rat_sub(&r->q, &a->q, &b->q); }
    r->m = m;
    normalize(r);
    return true;
}

bool qs_add(QSurd* r, const QSurd* a, const QSurd* b) { return add_signed(r, a, b, 1); }
bool qs_sub(QSurd* r, const QSurd* a, const QSurd* b) { return add_signed(r, a, b, -1); }

// (p1 + q1·√m)(p2 + q2·√m) = p1·p2 + q1·q2·m + (p1·q2 + q1·p2)·√m
bool qs_mul(QSurd* r, const QSurd* a, const QSurd* b) {
    long m = field_of(a, b);
    if (m < 0) return false;
    Rational p, q, t, mm = rat_from_int(m);
    rat_init(&p); rat_init(&q); rat_init(&t);
    rat_mul(&p, &a->p, &b->p);
    rat_mul(&t, &a->q, &b->q);
    rat_mul(&t, &t, &mm);
    rat_add(&p, &p, &t);
    rat_mul(&q, &a->p, &b->q);
    rat_mul(&t, &a->q, &b->p);
    rat_add(&q, &q, &t);
    rat_set(&r->p, &p); rat_set(&r->q, &q);
    r->m = m;
    normalize(r);
    rat_clear(&p); rat_clear(&q); rat_clear(&t); rat_clear(&mm);
    return true;
}

// a / b = a·conj(b) / N(b), N(b) = p² - q²·m (nonzero for b ≠ 0, as m is not a square)
bool qs_div(QSurd* r, const QSurd* a, const QSurd* b) {
    if (field_of(a, b) < 0) return false;
    Rational n, t, mm = rat_from_int(b->m);
    rat_init(&n); rat_init(&t);
    rat_mul(&n, &b->p, &b->p);
    rat_mul(&t, &b->q, &b->q);
    rat_mul(&t, &t, &mm);
    rat_sub(&n, &n, &t);
    bool ok = !rat_is_zero(&n);
    if (ok) {
        QSurd conj;
        qs_init(&conj);
        rat_set(&conj.p, &b->p);
        rat_neg(&conj.q, &b->q);
        conj.m = b->m;
        qs_mul(r, a, &conj);
        rat_div(&r->p, &r->p, &n);
        rat_div(&r->q, &r->q, &n);
        qs_clear(&conj);
    }
    rat_clear(&n); rat_clear(&t); rat_clear(&mm);
    return ok;
}

// sign of p + q·√m: the sign shared by p and q, else that of the larger of p², q²·m
int qs_sign(const QSurd* a) {
    int sp = mpq_sgn(a->p.val), sq = mpq_sgn(a->q.val);
    if (!sq || sp == sq) return sp ? sp : sq;
    if (!sp) return sq;
    Rational p2, q2m, mm = rat_from_int(a->m);
    rat_init(&p2); rat_init(&q2m);
    rat_mul(&p2, &a->p, &a->p);
    rat_mul(&q2m, &a->q, &a->q);
    rat_mul(&q2m, &q2m, &mm);
    int s = rat_cmp(&p2, &q2m) > 0 ? sp : sq;
    rat_clear(&p2); rat_clear(&q2m); rat_clear(&mm);
    return s;
}

static bool integral(const Rational* r) { return mpz_cmp_ui(mpq_denref(r->val), 1) == 0; }

// a + b·√m from integers
static void set_ints(QSurd* h, const mpz_t a, const mpz_t b, long m) {
    mpq_set_z(h->p.val, a);
    mpq_set_z(h->q.val, b);
    h->m = m;
    normalize(h);
}

// (a + b√m)² = p + q√m means a² + m·b² = p and 2ab = q, so with n² = p² - m·q² (the norm),
// a² = (p ± n)/2; each candidate is checked by squaring
bool qs_sqrt(QSurd* h, const QSurd* g) {
    if (!g->m || !integral(&g->p) || !integral(&g->q)) return false;
    mpz_t N, t, n, a, b;
    mpz_inits(N, t, n, a, b, NULL);
    mpz_srcptr p = mpq_numref(g->p.val), q = mpq_numref(g->q.val);
    mpz_mul(N, p, p);
    mpz_mul(t, q, q);
    mpz_mul_si(t, t, g->m);
    mpz_sub(N, N, t);
    bool found = false;
    if (mpz_sgn(N) >= 0 && mpz_perfect_square_p(N)) {
        mpz_sqrt(n, N);
        for (int s = 1; s >= -1 && !found; s -= 2) {
            if (s > 0) mpz_add(t, p, n); else mpz_sub(t, p, n);
            if (mpz_sgn(t) <= 0 || !mpz_even_p(t)) continue;
            mpz_divexact_ui(t, t, 2);
            if (!mpz_perfect_square_p(t)) continue;
            mpz_sqrt(a, t);
            mpz_mul_ui(t, a, 2);
            if (!mpz_divisible_p(q, t)) continue;
            mpz_divexact(b, q, t);
            QSurd sq;
            qs_init(&sq);
            set_ints(h, a, b, g->m);
            qs_mul(&sq, h, h);
            found = rat_eq(&sq.p, &g->p) && rat_eq(&sq.q, &g->q);
            qs_clear(&sq);
        }
        if (found && qs_sign(h) < 0) { rat_neg(&h->p, &h->p); rat_neg(&h->q, &h->q); }
    }
    mpz_clears(N, t, n, a, b, NULL);
    return found;
}

Expr* qs_expr(const QSurd* a) {
    Expr* surd = NULL;
    if (!rat_is_zero(&a->q)) {
        Expr* s = create_sqrt(create_const(a->m));
        surd = rat_is_one(&a->q) ? s : create_mul(create_const_rat(&a->q), s);
    }
    if (!surd) return create_const_rat(&a->p);
    return rat_is_zero(&a->p) ? surd : create_add(create_const_rat(&a->p), surd);
}
