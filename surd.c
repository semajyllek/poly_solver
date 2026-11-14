// Exact real constants: rationals times canonical radicals times trig values at rational
// multiples of π. Trig values on the table (multiples of π/6 and π/4) become rationals
// and square roots; the rest stay symbolic as cos/sin/cot/csc(kπ/N) with kπ/N in (0, π/2).
#include "surd.h"
#include <assert.h>
#include <stdlib.h>

#define TRIAL_LIMIT 10000UL

// --- radicals ---

// move the q-th powers found in N into s: afterwards N_old = s^q · N
static void pull_powers(mpz_t N, mpz_t s, unsigned long q) {
    mpz_t pq, r;
    mpz_inits(pq, r, NULL);
    mpz_set_ui(s, 1);
    for (unsigned long p = 2; p < TRIAL_LIMIT; p++) {
        mpz_ui_pow_ui(pq, p, q);
        if (mpz_cmp(pq, N) > 0) break;
        while (mpz_divisible_p(N, pq)) { mpz_divexact(N, N, pq); mpz_mul_ui(s, s, p); }
    }
    if (mpz_root(r, N, q)) { mpz_mul(s, s, r); mpz_set_ui(N, 1); } // a large power trial division missed
    mpz_clears(pq, r, NULL);
}

// N = m^j with j as large as possible (m = N, j = 1 when N is not a perfect power)
static unsigned long perfect_base(mpz_t m, const mpz_t N) {
    if (mpz_cmp_ui(N, 1) > 0)
        for (unsigned long k = mpz_sizeinbase(N, 2); k >= 2; k--)
            if (mpz_root(m, N, k)) return k;
    mpz_set(m, N);
    return 1;
}

static void rat_mul_z(Rational* r, const mpz_t z) {
    Rational t;
    rat_init(&t);
    mpq_set_z(t.val, z);
    rat_mul(r, r, &t);
    rat_clear(&t);
}

// x's radical becomes N^(1/v) = m^(j/v), with any whole power moved into c
static void set_radical(Mono* x, const mpz_t N, unsigned long v) {
    unsigned long j = perfect_base(x->m, N);
    unsigned long w = j / v;                              // nonzero only if trial division missed a power
    if (w > 0) {
        mpz_t mw;
        mpz_init(mw);
        mpz_pow_ui(mw, x->m, w);
        rat_mul_z(&x->c, mw);
        mpz_clear(mw);
    }
    rat_set_si(&x->e, (long)(j - w * v), (long)v);
    if (rat_is_zero(&x->e) || mpz_cmp_ui(x->m, 1) == 0) { mpz_set_ui(x->m, 1); rat_set_si(&x->e, 0, 1); }
}

void mono_init_radical(Mono* x, const Rational* c, const Power* pw) {
    const Rational *R = pw->R, *e = pw->e;
    rat_init_set(&x->c, c);
    mpz_init_set_ui(x->m, 1);
    rat_init(&x->e);
    rat_init(&x->angle);
    x->s = 1;
    x->kind = TRIG_NONE;
    mpz_t w, rem, N, b, t, s;
    mpz_inits(w, rem, N, b, t, s, NULL);
    mpz_fdiv_qr(w, rem, mpq_numref(e->val), mpq_denref(e->val)); // e = w + rem/v, 0 <= rem < v
    Rational Rw;
    rat_init(&Rw);
    rat_pow_int(&Rw, R, (int)mpz_get_si(w));
    rat_mul(&x->c, &x->c, &Rw);
    unsigned long r = mpz_get_ui(rem), v = mpz_get_ui(mpq_denref(e->val));
    mpz_pow_ui(N, mpq_numref(R->val), r);                  // (a/b)^(r/v) = (a^r · b^(r(v-1)))^(1/v) / b^r
    mpz_pow_ui(b, mpq_denref(R->val), r);
    mpz_pow_ui(t, b, v - 1);
    mpz_mul(N, N, t);
    pull_powers(N, s, v);
    rat_mul_z(&x->c, s);
    mpq_set_z(Rw.val, b);
    rat_div(&x->c, &x->c, &Rw);
    set_radical(x, N, v);
    rat_clear(&Rw);
    mpz_clears(w, rem, N, b, t, s, NULL);
}

void mono_clear(Mono* x) {
    rat_clear(&x->c);
    mpz_clear(x->m);
    rat_clear(&x->e);
    rat_clear(&x->angle);
}

// --- trig values ---

// a <- a mod period, in [0, period)
static void reduce_mod(Rational* a, long period) {
    Rational p, q;
    rat_init(&p); rat_init(&q);
    rat_set_si(&p, period, 1);
    rat_div(&q, a, &p);
    mpz_fdiv_q(mpq_numref(q.val), mpq_numref(q.val), mpq_denref(q.val)); // floor(a / period)
    mpz_set_ui(mpq_denref(q.val), 1);
    rat_mul(&q, &q, &p);
    rat_sub(a, a, &q);
    rat_clear(&p); rat_clear(&q);
}

// a > b/den
static bool above(const Rational* a, long b, long den) {
    Rational t = rat_from_si(b, den);
    bool r = rat_cmp(a, &t) > 0;
    rat_clear(&t);
    return r;
}

// a <- k - a
static void reflect(Rational* a, long k) {
    Rational t = rat_from_int(k);
    rat_sub(a, &t, a);
    rat_clear(&t);
}

// trig(a·π) = sign · trig(a'·π) with a' in [0, 1/2]
static int fold(TrigKind kind, Rational* a) {
    int sign = 1;
    if (kind == TRIG_COT) {
        reduce_mod(a, 1);
        if (above(a, 1, 2)) { sign = -1; reflect(a, 1); }
        return sign;
    }
    reduce_mod(a, 2);
    if (kind == TRIG_COS) {
        if (above(a, 1, 1)) reflect(a, 2);
        if (above(a, 1, 2)) { sign = -1; reflect(a, 1); }
    } else {                                   // sin, csc
        Rational one = rat_one_val();
        if (rat_cmp(a, &one) >= 0) { sign = -1; rat_sub(a, a, &one); }
        rat_clear(&one);
        if (above(a, 1, 2)) reflect(a, 1);
    }
    return sign;
}

typedef struct { int num, den, s; } TrigEntry; // num/den · sqrt(s); den = 0: undefined

// angles 0, π/6, π/4, π/3, π/2
static const TrigEntry TABLE[5][5] = {
    [TRIG_COS] = { { 1, 1, 1 }, { 1, 2, 3 }, { 1, 2, 2 }, { 1, 2, 1 }, { 0, 1, 1 } },
    [TRIG_SIN] = { { 0, 1, 1 }, { 1, 2, 1 }, { 1, 2, 2 }, { 1, 2, 3 }, { 1, 1, 1 } },
    [TRIG_COT] = { { 0, 0, 0 }, { 1, 1, 3 }, { 1, 1, 1 }, { 1, 3, 3 }, { 0, 1, 1 } },
    [TRIG_CSC] = { { 0, 0, 0 }, { 2, 1, 1 }, { 1, 1, 2 }, { 2, 3, 3 }, { 1, 1, 1 } },
};

// table row for a folded angle, or -1 when the value stays symbolic
static int table_index(const Rational* a) {
    Rational t = rat_from_int(12);
    rat_mul(&t, &t, a);
    int idx = -1;
    if (mpz_cmp_ui(mpq_denref(t.val), 1) == 0) {
        long k = mpz_get_si(mpq_numref(t.val));
        idx = k == 0 ? 0 : k == 2 ? 1 : k == 3 ? 2 : k == 4 ? 3 : k == 6 ? 4 : -1;
    }
    rat_clear(&t);
    return idx;
}

// x <- x · sqrt(s), merging with the radical when it is the same base
static void mul_sqrt(Mono* x, int s) {
    assert(x->s == 1);
    if (mpz_cmp_ui(x->m, (unsigned long)s) != 0) { x->s = s; return; }
    Rational half = rat_from_si(1, 2);
    rat_add(&x->e, &x->e, &half);
    rat_clear(&half);
    if (above(&x->e, 1, 1) || rat_is_one(&x->e)) {      // m^(e + 1/2) with e + 1/2 >= 1
        rat_mul_z(&x->c, x->m);
        Rational one = rat_one_val();
        rat_sub(&x->e, &x->e, &one);
        rat_clear(&one);
    }
    if (rat_is_zero(&x->e)) mpz_set_ui(x->m, 1);
}

void mono_mul_trig(Mono* x, TrigKind kind, const Rational* angle) {
    Rational a;
    rat_init_set(&a, angle);
    Rational sign = rat_from_int(fold(kind, &a));
    rat_mul(&x->c, &x->c, &sign);
    int idx = table_index(&a);
    if (idx < 0) {                                        // off the table: keep it symbolic
        assert(x->kind == TRIG_NONE);
        x->kind = kind;
        rat_set(&x->angle, &a);
    } else {
        TrigEntry v = TABLE[kind][idx];
        assert(v.den != 0);                               // cot/csc of a multiple of π
        Rational q = rat_from_si(v.num, v.den);
        rat_mul(&x->c, &x->c, &q);
        rat_clear(&q);
        if (v.s > 1) mul_sqrt(x, v.s);
    }
    rat_clear(&a);
    rat_clear(&sign);
}

// --- sums of monomials ---

RealConst rc_zero(void) { return (RealConst){ NULL, 0 }; }

bool rc_is_zero(const RealConst* k) { return k->count == 0; }

static bool like(const Mono* a, const Mono* b) {
    return mpz_cmp(a->m, b->m) == 0 && rat_eq(&a->e, &b->e) && a->s == b->s && a->kind == b->kind
        && (a->kind == TRIG_NONE || rat_eq(&a->angle, &b->angle));
}

void rc_add(RealConst* k, Mono* x) {
    for (int i = 0; i < k->count && !rat_is_zero(&x->c); i++) {
        if (!like(&k->t[i], x)) continue;
        rat_add(&k->t[i].c, &k->t[i].c, &x->c);
        mono_clear(x);
        if (rat_is_zero(&k->t[i].c)) { mono_clear(&k->t[i]); k->t[i] = k->t[--k->count]; }
        return;
    }
    if (rat_is_zero(&x->c)) { mono_clear(x); return; }
    k->t = realloc(k->t, (size_t)(k->count + 1) * sizeof(Mono));
    k->t[k->count++] = *x;                                // moves the GMP handles
}

void rc_clear(RealConst* k) {
    for (int i = 0; i < k->count; i++) mono_clear(&k->t[i]);
    free(k->t);
    k->t = NULL;
    k->count = 0;
}

// --- expressions ---

static Expr* mul_front(Expr* f, Expr* rest) { return rest ? create_mul(f, rest) : f; }

static Expr* add(Expr* acc, Expr* t) { return acc ? create_add(acc, t) : t; }

// m^e, printed as sqrt(m) when e = 1/2
static Expr* radical_expr(const Mono* x) {
    Rational m;
    rat_init(&m);
    mpq_set_z(m.val, x->m);
    Expr* base = create_const_rat(&m);
    rat_clear(&m);
    bool half = mpz_cmp_ui(mpq_denref(x->e.val), 2) == 0;
    return half ? create_sqrt(base) : create_pow(base, create_const_rat(&x->e));
}

// c · m^e · sqrt(s) · trig · body (body may be NULL)
static Expr* mono_expr(const Mono* x, Expr* body) {
    Expr* r = body;
    if (x->kind != TRIG_NONE) r = mul_front(create_trig(x->kind, &x->angle), r);
    if (x->s > 1) r = mul_front(create_sqrt(create_const(x->s)), r);
    if (mpz_cmp_ui(x->m, 1) != 0) r = mul_front(radical_expr(x), r);
    Expr* c = create_const_rat(&x->c);
    return r ? create_mul(c, r) : c;
}

Expr* rc_x_pow(int power) {
    Polynomial* p = create_polynomial();
    Rational one = rat_one_val();
    add_term(p, &one, power);
    rat_clear(&one);
    return create_poly_expr(p);
}

Expr* rc_add_terms(Expr* acc, const RealConst* k, int power) {
    for (int i = 0; i < k->count; i++) acc = add(acc, mono_expr(&k->t[i], power < 0 ? NULL : rc_x_pow(power)));
    return acc;
}

Expr* rc_scale(const RealConst* k, Expr* body) {
    if (k->count == 0) { free_expr(body); return NULL; }
    if (k->count == 1) return mono_expr(&k->t[0], body);
    return create_mul(rc_add_terms(NULL, k, -1), body);
}
