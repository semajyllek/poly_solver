// Exact folding of constant expressions produced by evaluating an antiderivative at a
// rational point. A term is c · f_1 · ... · f_k with c rational and f_i non-rational
// factors (radicals, sqrt, trig constants, ln, atan, π); a sum is a rational part, a
// multiple of π, logs of rationals c·ln(v) merged by v, and other terms merged when their
// factors print the same.
#include "constfold.h"
#include "surd.h"
#include "qsurd.h"
#include <stdlib.h>
#include <string.h>

static bool is_const(const Expr* e) { return e->type == EXPR_CONST; }

// --- products ---

typedef struct {
    Expr** f;
    int n;
} Factors;

static void push(Factors* F, Expr* e) {
    F->f = realloc(F->f, (size_t)(F->n + 1) * sizeof(Expr*));
    F->f[F->n++] = e;
}

static void free_factors(Factors* F) {
    for (int i = 0; i < F->n; i++) free_expr(F->f[i]);
    free(F->f);
    F->f = NULL;
    F->n = 0;
}

static bool is_sqrt_const(const Expr* e) { return e->type == EXPR_SQRT && is_const(e->fn_arg); }

// c · F = e, with c rational and F the non-rational factors (copies)
static void gather(const Expr* e, Rational* c, Factors* F) {
    if (is_const(e)) { rat_mul(c, c, &e->constant); return; }
    if (e->type == EXPR_NEG) { rat_neg(c, c); gather(e->operand, c, F); return; }
    if (e->type == EXPR_BINOP && e->binop.op == OP_MUL) {
        gather(e->binop.left, c, F);
        gather(e->binop.right, c, F);
        return;
    }
    const Expr* d = e->type == EXPR_BINOP && e->binop.op == OP_DIV ? e->binop.right : NULL;
    if (d && is_const(d) && !rat_is_zero(&d->constant)) { // x / q = (1/q)·x
        rat_div(c, c, &d->constant);
        gather(e->binop.left, c, F);
        return;
    }
    if (d && is_sqrt_const(d)) {                         // x / sqrt(m) = (1/m)·x·sqrt(m)
        rat_div(c, c, &d->fn_arg->constant);
        gather(e->binop.left, c, F);
        push(F, expr_copy(d));
        return;
    }
    push(F, expr_copy(e));
}

// the square roots of positive integers in F multiply into one, squares moved into c:
// sqrt(m)·sqrt(m) = m, sqrt(2)·sqrt(3) = sqrt(6), sqrt(2)·sqrt(6) = 2·sqrt(3)
static bool is_int_sqrt(const Expr* e) {
    return e && is_sqrt_const(e) && mpz_cmp_ui(mpq_denref(e->fn_arg->constant.val), 1) == 0
        && mpz_sgn(mpq_numref(e->fn_arg->constant.val)) > 0;
}

static void pair_sqrts(Rational* c, Factors* F) {
    mpz_t N, q, r;
    mpz_inits(N, q, r, NULL);
    mpz_set_ui(N, 1);
    int n = 0;
    for (int i = 0; i < F->n; i++)
        if (is_int_sqrt(F->f[i])) { mpz_mul(N, N, mpq_numref(F->f[i]->fn_arg->constant.val)); free_expr(F->f[i]); F->f[i] = NULL; n++; }
    mpz_set_ui(q, 1);                                     // N = q²·r, r squarefree by trial division
    for (unsigned long p = 2; p * p <= 1000000 && mpz_cmp_ui(N, p * p) >= 0; p++)
        while (mpz_divisible_ui_p(N, p * p)) { mpz_divexact_ui(N, N, p * p); mpz_mul_ui(q, q, p); }
    if (mpz_perfect_square_p(N)) { mpz_sqrt(r, N); mpz_mul(q, q, r); mpz_set_ui(N, 1); }
    Rational t;
    rat_init(&t);
    mpq_set_z(t.val, q);
    rat_mul(c, c, &t);
    int k = 0;
    for (int i = 0; i < F->n; i++) if (F->f[i]) F->f[k++] = F->f[i];
    F->n = k;
    if (n && mpz_cmp_ui(N, 1) != 0) {                     // the remaining sqrt(r)
        mpq_set_z(t.val, N);
        F->f[F->n++] = create_sqrt(create_const_rat(&t));
    }
    rat_clear(&t);
    mpz_clears(N, q, r, NULL);
}

// print order of a factor: radicals, trig values, π, then ln/atan and anything else
static int factor_rank(const Expr* e) {
    if (e->type == EXPR_SQRT || (e->type == EXPR_BINOP && e->binop.op == OP_POW)) return 0;
    if (e->type == EXPR_TRIG) return e->trig.kind == TRIG_NONE ? 2 : 1;
    return 3;
}

// stable insertion sort of the factors by rank
static void order_factors(Factors* F) {
    for (int i = 1; i < F->n; i++)
        for (int j = i; j > 0 && factor_rank(F->f[j - 1]) > factor_rank(F->f[j]); j--) {
            Expr* t = F->f[j]; F->f[j] = F->f[j - 1]; F->f[j - 1] = t;
        }
}

// c · Π F as an expression, factors in print order; takes F's factors
static Expr* make_product(const Rational* c, Factors* F) {
    if (rat_is_zero(c)) { free_factors(F); return create_const(0); }
    order_factors(F);
    Expr* body = NULL;
    for (int i = F->n - 1; i >= 0; i--) body = body ? create_mul(F->f[i], body) : F->f[i];
    free(F->f);
    F->f = NULL;
    F->n = 0;
    if (!body) return create_const_rat(c);
    return rat_is_one(c) ? body : create_mul(create_const_rat(c), body);
}

// l · r folded; takes l and r
static Expr* fold_product(Expr* l, Expr* r) {
    Rational c = rat_one_val();
    Factors F = { NULL, 0 };
    gather(l, &c, &F);
    gather(r, &c, &F);
    free_expr(l); free_expr(r);
    pair_sqrts(&c, &F);
    Expr* e = make_product(&c, &F);
    rat_clear(&c);
    return e;
}

// --- sums ---

typedef struct {
    Rational c;
    Expr* body; // product of non-rational factors
    char* key;  // its printed form
} SumTerm;

typedef struct {
    Rational rational, pi;
    Rational *lc, *lv; // lc[i]·ln(lv[i]), lv[i] > 1
    int nl;
    SumTerm* t;
    int nt;
} Sum;

// the smallest k >= 2 for which v > 0 is a perfect k-th power, with n/d set to the roots; or 0
static long root_degree(const Rational* v, mpz_t n, mpz_t d) {
    size_t nb = mpz_sizeinbase(mpq_numref(v->val), 2), db = mpz_sizeinbase(mpq_denref(v->val), 2);
    long top = (long)(nb > db ? nb : db);                 // 2^k <= v needs k <= bits
    for (long k = 2; k <= top; k++)
        if (mpz_root(n, mpq_numref(v->val), (unsigned long)k) && mpz_root(d, mpq_denref(v->val), (unsigned long)k)) return k;
    return 0;
}

// v = w^k with k as large as possible, for v > 0: v becomes w, and k is returned
static long perfect_root(Rational* v) {
    long total = 1, k;
    mpz_t n, d;
    mpz_init(n); mpz_init(d);
    while ((k = root_degree(v, n, d))) {
        mpq_set_num(v->val, n);
        mpq_set_den(v->val, d);
        total *= k;
    }
    mpz_clear(n); mpz_clear(d);
    return total;
}

static Sum sum_new(void) { return (Sum){ rat_zero_val(), rat_zero_val(), NULL, NULL, 0, NULL, 0 }; }

// c·ln(v), merged with the group of the same v (written with v > 1)
static void add_ln(Sum* S, Rational* c, Rational* v) {
    Rational one = rat_one_val();
    if (rat_cmp(v, &one) < 0) { mpq_inv(v->val, v->val); rat_neg(c, c); }
    bool unit = rat_eq(v, &one);
    rat_clear(&one);
    if (unit) return;
    for (int i = 0; i < S->nl; i++)
        if (rat_eq(&S->lv[i], v)) { rat_add(&S->lc[i], &S->lc[i], c); return; }
    S->lc = realloc(S->lc, (size_t)(S->nl + 1) * sizeof(Rational));
    S->lv = realloc(S->lv, (size_t)(S->nl + 1) * sizeof(Rational));
    rat_init_set(&S->lc[S->nl], c);
    rat_init_set(&S->lv[S->nl], v);
    S->nl++;
}

// c·body, merged with a like term
static void add_other(Sum* S, const Rational* c, Expr* body) {
    char* key = expr_to_string(body);
    for (int i = 0; i < S->nt; i++)
        if (strcmp(S->t[i].key, key) == 0) {
            rat_add(&S->t[i].c, &S->t[i].c, c);
            free(key);
            free_expr(body);
            return;
        }
    S->t = realloc(S->t, (size_t)(S->nt + 1) * sizeof(SumTerm));
    rat_init_set(&S->t[S->nt].c, c);
    S->t[S->nt].body = body;
    S->t[S->nt].key = key;
    S->nt++;
}

static bool is_sum(const Expr* e) {
    return e->type == EXPR_BINOP && (e->binop.op == OP_ADD || e->binop.op == OP_SUB);
}

static void add_summand(Sum* S, const Expr* e);

// coef · (each summand of s), added one by one
static void distribute(Sum* S, const Expr* coef, const Expr* s) {
    if (is_sum(s)) {
        distribute(S, coef, s->binop.left);
        if (s->binop.op == OP_ADD) { distribute(S, coef, s->binop.right); return; }
        Expr* neg = create_mul(create_const(-1), expr_copy(coef));
        distribute(S, neg, s->binop.right);
        free_expr(neg);
        return;
    }
    Expr* p = create_mul(expr_copy(coef), expr_copy(s));
    add_summand(S, p);
    free_expr(p);
}

// a product containing a sum factor is expanded over it; true if it was
static bool expand_product(Sum* S, Rational* c, Factors* F) {
    for (int i = 0; i < F->n; i++) {
        if (!is_sum(F->f[i])) continue;
        Expr* s = F->f[i];
        F->f[i] = create_const(1);
        Expr* coef = make_product(c, F);
        distribute(S, coef, s);
        free_expr(coef);
        free_expr(s);
        return true;
    }
    return false;
}

// one summand (not a sum itself), classified
static void add_summand(Sum* S, const Expr* e) {
    Rational c = rat_one_val();
    Factors F = { NULL, 0 };
    gather(e, &c, &F);
    pair_sqrts(&c, &F);
    if (expand_product(S, &c, &F)) { rat_clear(&c); return; }
    const Expr* f = F.n == 1 ? F.f[0] : NULL;
    if (F.n == 0) rat_add(&S->rational, &S->rational, &c);
    else if (f && f->type == EXPR_TRIG && f->trig.kind == TRIG_NONE) rat_add(&S->pi, &S->pi, &c);
    else if (f && f->type == EXPR_LN && is_const(f->ln_arg)) {
        Rational v;
        rat_init_set(&v, &f->ln_arg->constant);
        add_ln(S, &c, &v);
        rat_clear(&v);
    } else {
        Rational one = rat_one_val();
        add_other(S, &c, make_product(&one, &F));
        rat_clear(&one);
    }
    free_factors(&F);
    rat_clear(&c);
}

// every summand of e, with sign applied
static void add_all(Sum* S, const Expr* e, int sign) {
    bool add = e->type == EXPR_BINOP && e->binop.op == OP_ADD, sub = e->type == EXPR_BINOP && e->binop.op == OP_SUB;
    if (add || sub) {
        add_all(S, e->binop.left, sign);
        add_all(S, e->binop.right, sub ? -sign : sign);
    } else if (e->type == EXPR_NEG) {
        add_all(S, e->operand, -sign);
    } else if (sign > 0) {
        add_summand(S, e);
    } else {
        Expr* n = create_mul(create_const(-1), expr_copy(e));
        add_summand(S, n);
        free_expr(n);
    }
}

static Expr* append(Expr* acc, Expr* t) { return acc ? create_add(acc, t) : t; }

// c·body; takes body
static Expr* scaled(const Rational* c, Expr* body) {
    Factors F = { NULL, 0 };
    push(&F, body);
    return make_product(c, &F);
}

// c_j / c_i as a small integer k (so c_i·ln(a) + c_j·ln(b) = c_i·ln(a·b^k)), or 0
static long log_ratio(const Rational* ci, const Rational* cj) {
    if (rat_is_zero(ci) || rat_is_zero(cj)) return 0;
    Rational q;
    rat_init(&q);
    rat_div(&q, cj, ci);
    bool small = mpz_cmp_ui(mpq_denref(q.val), 1) == 0 && mpz_cmpabs_ui(mpq_numref(q.val), 64) <= 0;
    long k = small ? mpz_get_si(mpq_numref(q.val)) : 0;
    rat_clear(&q);
    return k;
}

// fold log j into log i: v_i <- v_i · v_j^k, c_j <- 0
static void merge_log(Sum* S, int i, int j) {
    long k = log_ratio(&S->lc[i], &S->lc[j]);
    Rational p;
    rat_init(&p);
    rat_pow_int(&p, &S->lv[j], (int)k);
    rat_mul(&S->lv[i], &S->lv[i], &p);
    rat_set_si(&S->lc[j], 0, 1);
    rat_clear(&p);
}

// merge logs whose coefficients are integer multiples of one another, then write v > 1
static void merge_logs(Sum* S) {
    for (int i = 0; i < S->nl; i++)
        for (int j = 0; j < S->nl; j++) {
            if (i == j || rat_is_zero(&S->lc[i])) continue;
            if (log_ratio(&S->lc[i], &S->lc[j])) merge_log(S, i, j);
        }
    Rational one = rat_one_val();
    for (int i = 0; i < S->nl; i++)
        if (rat_cmp(&S->lv[i], &one) < 0) { mpq_inv(S->lv[i].val, S->lv[i].val); rat_neg(&S->lc[i], &S->lc[i]); }
        else if (rat_eq(&S->lv[i], &one)) rat_set_si(&S->lc[i], 0, 1);
    rat_clear(&one);
}

// --- ln and atan terms over Q(sqrt(m)) ---

static Expr* fold_ln(Expr* a);
static Expr* fold_atan(Expr* a);

// a summand body rest·fn(arg) with fn = ln or atan and arg in Q(sqrt(m))
typedef struct {
    ExprType type; // EXPR_LN or EXPR_ATAN
    Expr* rest;    // the other factors
    char* key;     // rest, printed
    QSurd arg;
} FnTerm;

static bool fn_term(const Expr* body, FnTerm* t) {
    Rational c = rat_one_val();
    Factors F = { NULL, 0 };
    gather(body, &c, &F);
    int k = -1, fns = 0;
    for (int i = 0; i < F.n; i++)
        if (F.f[i]->type == EXPR_LN || F.f[i]->type == EXPR_ATAN) { k = i; fns++; }
    qs_init(&t->arg);
    bool ok = fns == 1 && qs_parse(F.f[k]->type == EXPR_LN ? F.f[k]->ln_arg : F.f[k]->fn_arg, &t->arg);
    if (ok) {
        t->type = F.f[k]->type;
        free_expr(F.f[k]);
        F.f[k] = F.f[--F.n];                              // drop the fn factor
        t->rest = make_product(&c, &F);
        t->key = expr_to_string(t->rest);
    } else {
        free_factors(&F);
        qs_clear(&t->arg);
    }
    rat_clear(&c);
    return ok;
}

static void fn_term_free(FnTerm* t) { free_expr(t->rest); free(t->key); qs_clear(&t->arg); }

// two like terms c·R·fn(u) and (p/q)·c·R·fn(v), merged as (c/q)·R·(q·fn(u) + p·fn(v))
typedef struct {
    const FnTerm *a, *b;
    long p, q; // small: |p|, q <= 64
} FnPair;

// c_j / c_i = p/q in lowest terms with |p|, q <= 64; false otherwise
static bool small_ratio(const Rational* ci, const Rational* cj, long* p, long* q) {
    if (rat_is_zero(ci) || rat_is_zero(cj)) return false;
    Rational r;
    rat_init(&r);
    rat_div(&r, cj, ci);
    bool small = mpz_cmpabs_ui(mpq_numref(r.val), 64) <= 0 && mpz_cmp_ui(mpq_denref(r.val), 64) <= 0;
    if (small) { *p = mpz_get_si(mpq_numref(r.val)); *q = mpz_get_si(mpq_denref(r.val)); }
    rat_clear(&r);
    return small;
}

// g ← g·a^n (n may be negative)
static bool mul_pow(QSurd* g, const QSurd* a, long n) {
    bool ok = true;
    for (long i = 0; ok && i < labs(n); i++) ok = n > 0 ? qs_mul(g, g, a) : qs_div(g, g, a);
    return ok;
}

// the positive rational r with g = r·g', g' having coprime integer parts
static void content(Rational* r, const QSurd* g) {
    mpz_t n, d;
    mpz_init(n); mpz_init(d);
    mpz_gcd(n, mpq_numref(g->p.val), mpq_numref(g->q.val));
    mpz_lcm(d, mpq_denref(g->p.val), mpq_denref(g->q.val));
    mpq_set_num(r->val, n);
    mpq_set_den(r->val, d);
    mpq_canonicalize(r->val);
    mpz_clear(n); mpz_clear(d);
}

// ln|g| = acc + sign·power·ln(g), with g reduced in place by canon_ln
typedef struct {
    Expr* acc; // rational logs split off so far
    long sign, power;
} LnForm;

// acc += sign·power·k·ln(v)
// acc += sign·power·k·ln|v|, written with |v| > 1 (ln(1/17) = -ln(17); ln(1) dropped)
static void split_off(LnForm* f, const Rational* k, const Rational* v) {
    Rational c = rat_from_int(f->sign * f->power), a, one = rat_one_val();
    rat_init(&a);
    rat_mul(&c, &c, k);
    rat_abs(&a, v);
    if (rat_cmp(&a, &one) < 0) { mpq_inv(a.val, a.val); rat_neg(&c, &c); }
    Rational power = rat_from_int(perfect_root(&a));      // ln(697²) = 2·ln(697)
    rat_mul(&c, &c, &power);
    rat_clear(&power);
    if (!rat_eq(&a, &one)) {
        Expr* t = create_mul(create_const_rat(&c), create_ln(create_const_rat(&a)));
        f->acc = f->acc ? create_add(f->acc, t) : t;
    }
    rat_clear(&c); rat_clear(&a); rat_clear(&one);
}

static bool integral(const Rational* r) { return mpz_cmp_ui(mpq_denref(r->val), 1) == 0; }

static bool below_one(const QSurd* g) {
    QSurd one, d;
    qs_init(&one); qs_init(&d);
    rat_set_si(&one.p, 1, 1);
    qs_sub(&d, g, &one);
    bool below = qs_sign(&d) < 0;
    qs_clear(&one); qs_clear(&d);
    return below;
}

// g/√m = q + (p/m)·√m when m divides the integer p: ln|g| = 1/2·ln(m) + ln|g/√m|
static bool split_sqrt_m(QSurd* g, LnForm* f) {
    bool divides = integral(&g->p) && mpz_divisible_ui_p(mpq_numref(g->p.val), (unsigned long)g->m);
    if (!divides) return false;
    Rational half = rat_from_si(1, 2), m = rat_from_int(g->m), p;
    rat_init_set(&p, &g->p);
    split_off(f, &half, &m);
    rat_set(&g->p, &g->q);
    rat_div(&g->q, &p, &m);
    if (rat_is_zero(&g->q)) g->m = 0;
    rat_clear(&half); rat_clear(&m); rat_clear(&p);
    return true;
}

// with g > 1 but p or q negative, ln g = ln|N| - ln|conj(g)|, N = p² - q²·m, and
// |conj(g)| = |p| + |q|·√m: conjugate arguments (1 + 3√2, -1 + 3√2) print alike
static bool split_conjugate(QSurd* g, LnForm* f) {
    if (mpq_sgn(g->p.val) > 0 && mpq_sgn(g->q.val) > 0) return false;
    Rational n, t, one = rat_one_val(), mm = rat_from_int(g->m);
    rat_init(&n); rat_init(&t);
    rat_mul(&n, &g->p, &g->p);
    rat_mul(&t, &g->q, &g->q);
    rat_mul(&t, &t, &mm);
    rat_sub(&n, &n, &t);
    split_off(f, &one, &n);
    f->sign = -f->sign;
    rat_abs(&g->p, &g->p);
    rat_abs(&g->q, &g->q);
    rat_clear(&n); rat_clear(&t); rat_clear(&one); rat_clear(&mm);
    return true;
}

// reduce a surd g to the canonical log argument: positive, > 1, content 1, no factor √m
// (each split lowers |norm|, so the loop ends), both parts positive (conjugate swap,
// at most once per norm), then square roots taken while g is a
// perfect square: ln|2 - sqrt(2)| = 1/2·ln(2) - ln(1 + sqrt(2)), ln(3 + 2·sqrt(2)) =
// 2·ln(1 + sqrt(2)). Returns the number of rewrites (0: g was canonical already).
static int canon_ln(QSurd* g, LnForm* f) {
    int steps = 0;
    Rational r, one = rat_one_val();
    rat_init(&r);
    while (g->m) {
        if (qs_sign(g) < 0) { rat_neg(&g->p, &g->p); rat_neg(&g->q, &g->q); steps++; }
        if (below_one(g)) { QSurd u; qs_init(&u); rat_set_si(&u.p, 1, 1); qs_div(g, &u, g); qs_clear(&u); f->sign = -f->sign; steps++; }
        content(&r, g);
        if (!rat_is_one(&r)) { split_off(f, &one, &r); rat_div(&g->p, &g->p, &r); rat_div(&g->q, &g->q, &r); steps++; }
        if (split_sqrt_m(g, f)) { steps++; continue; }
        if (split_conjugate(g, f)) { steps++; continue; }
        QSurd h;
        qs_init(&h);
        bool square = qs_sqrt(&h, g);
        if (square) { qs_set(g, &h); f->power *= 2; steps++; }
        qs_clear(&h);
        if (!square) break;
    }
    if (!g->m) { rat_abs(&g->p, &g->p); split_off(f, &one, &g->p); rat_set_si(&g->p, 1, 1); steps++; }
    rat_clear(&r); rat_clear(&one);
    return steps;
}

// g is already the canonical argument of its log
static bool ln_normal(const QSurd* g) {
    QSurd c;
    qs_init(&c);
    qs_set(&c, g);
    LnForm f = { NULL, 1, 1 };
    bool normal = canon_ln(&c, &f) == 0;
    free_expr(f.acc);
    qs_clear(&c);
    return normal;
}

// ln|g| for g in Q(sqrt(m)): the split-off rational logs plus sign·power·ln(canonical g)
static Expr* ln_surd(QSurd* g) {
    LnForm f = { NULL, 1, 1 };
    if (!g->m) {                                          // a rational: ln|v| with |v| > 1
        Rational one = rat_one_val();
        split_off(&f, &one, &g->p);
        rat_clear(&one);
        return f.acc ? f.acc : create_const(0);
    }
    canon_ln(g, &f);
    Expr* main = NULL;
    if (g->m) main = create_mul(create_const(f.sign * f.power), create_ln(qs_expr(g)));
    if (!f.acc) return main ? main : create_const(0);
    return main ? create_add(f.acc, main) : f.acc;
}

// q·ln|u| + p·ln|v| = ln|u^q·v^p|
static Expr* ln_merge(const FnPair* pr) {
    QSurd g;
    qs_init(&g);
    rat_set_si(&g.p, 1, 1);
    bool ok = mul_pow(&g, &pr->a->arg, pr->q) && mul_pow(&g, &pr->b->arg, pr->p);
    Expr* e = ok && qs_sign(&g) ? ln_surd(&g) : NULL;
    qs_clear(&g);
    return e;
}

// an angle atan(w) + n·π/2, built up exactly
typedef struct {
    QSurd w;
    long n;
} Angle;

// angle += atan(v): atan(w) + atan(v) = atan((w + v)/(1 - wv)), plus sign(w)·π when
// 1 - wv < 0; exactly sign(w)·π/2 (w = 0, n += sign(w)) when 1 - wv = 0
static bool add_atan(Angle* A, const QSurd* v) {
    QSurd wv, d, one;
    qs_init(&wv); qs_init(&d); qs_init(&one);
    rat_set_si(&one.p, 1, 1);
    int sw = qs_sign(&A->w);
    bool ok = qs_mul(&wv, &A->w, v) && qs_sub(&d, &one, &wv) && qs_add(&A->w, &A->w, v);
    int sd = ok ? qs_sign(&d) : 1;
    if (ok && sd == 0) { rat_set_si(&A->w.p, 0, 1); rat_set_si(&A->w.q, 0, 1); A->w.m = 0; A->n += sw; }
    else if (ok) { ok = qs_div(&A->w, &A->w, &d); if (sd < 0) A->n += 2 * sw; }
    qs_clear(&wv); qs_clear(&d); qs_clear(&one);
    return ok;
}

// q·atan(u) + p·atan(v) by repeated addition (|p| + q <= 128 steps), as atan(w) + n·π/2
static Expr* atan_merge(const FnPair* pr) {
    Angle A = { .n = 0 };
    QSurd v;
    qs_init(&A.w); qs_init(&v);
    qs_set(&v, &pr->b->arg);
    if (pr->p < 0) { rat_neg(&v.p, &v.p); rat_neg(&v.q, &v.q); }   // p·atan(v) = |p|·atan(-v)
    bool ok = true;
    for (long i = 0; ok && i < pr->q; i++) ok = add_atan(&A, &pr->a->arg);
    for (long i = 0; ok && i < labs(pr->p); i++) ok = add_atan(&A, &v);
    Expr* e = NULL;
    if (ok) {
        Rational turns = rat_from_si(A.n, 2);
        e = create_add(fold_atan(qs_expr(&A.w)), create_mul(create_const_rat(&turns), create_pi()));
        rat_clear(&turns);
    }
    qs_clear(&A.w); qs_clear(&v);
    return e;
}

// terms i and j folded into one when they are like ln or atan terms; true if they were
static bool try_merge(Sum* S, int i, int j) {
    if (rat_is_zero(&S->t[i].c) || rat_is_zero(&S->t[j].c)) return false;
    FnTerm a, b;
    if (!fn_term(S->t[i].body, &a)) return false;
    bool both = fn_term(S->t[j].body, &b);
    Expr* part = NULL;
    FnPair pr = { &a, &b, 0, 1 };
    if (both && a.type == b.type && strcmp(a.key, b.key) == 0 && small_ratio(&S->t[i].c, &S->t[j].c, &pr.p, &pr.q)) {
        bool same_kind = !a.arg.m == !b.arg.m;          // a merged surd log splits off rational logs,
        bool exact_ln = !a.arg.m || pr.q == 1;            // so mixing kinds would cycle; a surd log
        if (a.type == EXPR_LN && same_kind && exact_ln) part = ln_merge(&pr); // needs an integer ratio
        else if (a.type == EXPR_ATAN) part = atan_merge(&pr);
    }
    if (part) {
        Rational c = rat_from_si(1, pr.q);                // (c/q)·R·(q·fn(u) + p·fn(v))
        rat_mul(&c, &c, &S->t[i].c);
        Expr* e = scaled(&c, create_mul(expr_copy(a.rest), part));
        rat_clear(&c);
        rat_set_si(&S->t[i].c, 0, 1);
        rat_set_si(&S->t[j].c, 0, 1);
        add_all(S, e, 1);
        free_expr(e);
    }
    fn_term_free(&a);
    if (both) fn_term_free(&b);
    return part != NULL;
}

static bool is_power(const Rational* v) {
    mpz_t n, d;
    mpz_init(n); mpz_init(d);
    bool p = root_degree(v, n, d) > 0;
    mpz_clear(n); mpz_clear(d);
    return p;
}

// term i rewritten when it is a log not in canonical form: a surd as canon_ln makes it
// (ln|1 - sqrt(2)| = -ln(1 + sqrt(2))), a rational as v > 1; true if it was
static bool try_normalize(Sum* S, int i) {
    if (rat_is_zero(&S->t[i].c)) return false;
    FnTerm a;
    if (!fn_term(S->t[i].body, &a)) return false;
    Rational one = rat_one_val();
    bool rational_ok = !a.arg.m && rat_cmp(&a.arg.p, &one) > 0 && !is_power(&a.arg.p); // ln(17), not ln(1/17), ln(289)
    rat_clear(&one);
    bool go = a.type == EXPR_LN && (a.arg.m ? !ln_normal(&a.arg) : !rational_ok);
    if (go) {
        Expr* e = scaled(&S->t[i].c, create_mul(expr_copy(a.rest), ln_surd(&a.arg)));
        rat_set_si(&S->t[i].c, 0, 1);
        add_all(S, e, 1);
        free_expr(e);
    }
    fn_term_free(&a);
    return go;
}

// normalize logs, then merge like ln and atan terms pairwise until none combine. This
// ends: nothing creates an unnormalized log or a new surd log, each surd-log merge removes
// one surd log (it may add rational logs), and rational-log and atan merges each remove a
// term of their own kind without adding any other.
static void combine_fns(Sum* S) {
    for (int i = 0; i < S->nt; i++) try_normalize(S, i);  // appended terms are normal already
    for (bool again = true; again;) {
        again = false;
        for (int i = 0; i < S->nt && !again; i++)
            for (int j = 0; j < S->nt && !again; j++)     // both orders: surd logs need
                again = j != i && try_merge(S, i, j);     // c_j/c_i to be an integer
    }
}

// a product has a factor of this type (TRIG_NONE alone for π, when type is EXPR_TRIG)
static bool has_factor(const Expr* e, ExprType type) {
    if (e->type == type) return type != EXPR_TRIG || e->trig.kind == TRIG_NONE;
    return e->type == EXPR_BINOP && e->binop.op == OP_MUL && (has_factor(e->binop.left, type) || has_factor(e->binop.right, type));
}

// print order of a sum's terms: logs, atans, the rest, then multiples of π
static int term_rank(const Expr* body) {
    if (has_factor(body, EXPR_TRIG)) return 3;
    if (has_factor(body, EXPR_LN)) return 0;
    return has_factor(body, EXPR_ATAN) ? 1 : 2;
}

// by rank, then by printed form, so equal sums print alike whatever their build order
static int term_cmp(const void* x, const void* y) {
    const SumTerm *a = *(SumTerm* const*)x, *b = *(SumTerm* const*)y;
    int ra = term_rank(a->body), rb = term_rank(b->body);
    return ra != rb ? ra - rb : strcmp(a->key, b->key);
}

// the folded sum; the rational part goes first when rational_first, else last
static Expr* sum_expr(Sum* S, bool rational_first) {
    combine_fns(S);
    merge_logs(S);
    Expr* e = NULL;
    bool rat = !rat_is_zero(&S->rational);
    if (rat && rational_first) e = create_const_rat(&S->rational);
    SumTerm** order = malloc((size_t)(S->nt + 1) * sizeof(SumTerm*));
    int n = 0;
    for (int i = 0; i < S->nt; i++) if (!rat_is_zero(&S->t[i].c)) order[n++] = &S->t[i];
    qsort(order, (size_t)n, sizeof(SumTerm*), term_cmp);   // the same terms print alike
    for (int i = 0; i < n; i++) {
        e = append(e, scaled(&order[i]->c, order[i]->body));
        order[i]->body = NULL;                            // now owned by e
    }
    free(order);
    for (int i = 0; i < S->nt; i++) {
        free_expr(S->t[i].body);
        free(S->t[i].key);
        rat_clear(&S->t[i].c);
    }
    for (int i = 0; i < S->nl; i++) {
        Rational k = rat_from_int(perfect_root(&S->lv[i]));   // ln(4) = 2·ln(2)
        rat_mul(&S->lc[i], &S->lc[i], &k);
        rat_clear(&k);
        if (!rat_is_zero(&S->lc[i])) e = append(e, scaled(&S->lc[i], create_ln(create_const_rat(&S->lv[i]))));
        rat_clear(&S->lc[i]);
        rat_clear(&S->lv[i]);
    }
    if (!rat_is_zero(&S->pi)) e = append(e, scaled(&S->pi, create_pi()));
    if (rat && !rational_first) e = append(e, create_const_rat(&S->rational));
    if (!e) e = create_const(0);
    free(S->t); free(S->lc); free(S->lv);
    rat_clear(&S->rational); rat_clear(&S->pi);
    return e;
}

// l ± r folded; takes l and r
static Expr* fold_sum2(Expr* l, Expr* r, int sign) {
    Sum S = sum_new();
    add_all(&S, l, 1);
    add_all(&S, r, sign);
    free_expr(l); free_expr(r);
    return sum_expr(&S, true);
}

// --- ln and atan of constants ---

// the first printed summand of a sum starts with a minus sign
static bool leading_negative(const Expr* e) {
    while (is_sum(e)) e = e->binop.left;
    if (is_const(e)) return rat_is_negative(&e->constant);
    return e->type == EXPR_BINOP && e->binop.op == OP_MUL && is_const(e->binop.left)
        && rat_is_negative(&e->binop.left->constant);
}

// ln|f| for one non-rational factor: m^e -> e·ln(m), sqrt(m) -> 1/2·ln(m), else ln|f|
// (a sum is written with a positive leading term, since |-u| = |u|)
static Expr* ln_factor(Expr* f) {
    if (is_sum(f) && leading_negative(f)) f = fold_sum2(create_const(0), f, -1);
    bool radical = f->type == EXPR_BINOP && f->binop.op == OP_POW && is_const(f->binop.left) && is_const(f->binop.right);
    if (!radical && !is_sqrt_const(f)) return create_ln(f);
    Rational e, m;
    if (radical) { rat_init_set(&e, &f->binop.right->constant); rat_init_set(&m, &f->binop.left->constant); }
    else { e = rat_from_si(1, 2); rat_init_set(&m, &f->fn_arg->constant); }
    free_expr(f);
    Expr* r = create_mul(create_const_rat(&e), create_ln(create_const_rat(&m)));
    rat_clear(&e); rat_clear(&m);
    return r;
}

// ln|a| = ln|c| + Σ ln|f_i| for a = c·Π f_i, and ln|n/d| = ln|n| - ln|d|; takes a
static Expr* fold_ln(Expr* a) {
    if (a->type == EXPR_BINOP && a->binop.op == OP_DIV && !is_const(a->binop.right)) {
        Expr *n = fold_ln(expr_copy(a->binop.left)), *d = fold_ln(expr_copy(a->binop.right));
        free_expr(a);
        return fold_sum2(n, d, -1);
    }
    Rational c = rat_one_val();
    Factors F = { NULL, 0 };
    gather(a, &c, &F);
    free_expr(a);
    pair_sqrts(&c, &F);
    rat_abs(&c, &c);
    Expr* e = create_ln(create_const_rat(&c));
    for (int i = 0; i < F.n; i++) e = fold_sum2(e, ln_factor(F.f[i]), 1);
    free(F.f);
    rat_clear(&c);
    return fold_sum2(e, create_const(0), 1);          // drops ln(1) and merges the logs
}

// atan(c·sqrt(3)) for c = ±1, ±1/3 is ±π/3, ±π/6; atan(±1) = ±π/4; atan(0) = 0
static bool atan_value(const Rational* c, const Factors* F, Rational* turns) {
    Rational three = rat_from_int(3), third = rat_from_si(1, 3), a;
    rat_init(&a);
    rat_abs(&a, c);
    bool sqrt3 = F->n == 1 && is_sqrt_const(F->f[0]) && rat_eq(&F->f[0]->fn_arg->constant, &three);
    bool found = true;
    if (F->n == 0 && rat_is_zero(c)) rat_set_si(turns, 0, 1);
    else if (F->n == 0 && rat_is_one(&a)) rat_set_si(turns, 1, 4);
    else if (sqrt3 && rat_is_one(&a)) rat_set_si(turns, 1, 3);
    else if (sqrt3 && rat_eq(&a, &third)) rat_set_si(turns, 1, 6);
    else found = false;
    if (found && rat_is_negative(c)) rat_neg(turns, turns);
    rat_clear(&three); rat_clear(&third); rat_clear(&a);
    return found;
}

// atan(±(sqrt(2) ∓ 1)) and atan(±(2 ∓ sqrt(3))): tan of π/8, 3π/8, π/12, 5π/12
static bool atan_surd_value(const Expr* a, Rational* turns) {
    QSurd u;
    qs_init(&u);
    bool found = qs_parse(a, &u) && u.m && (u.m == 2 || u.m == 3);
    int sign = found ? qs_sign(&u) : 0;
    if (sign < 0) { rat_neg(&u.p, &u.p); rat_neg(&u.q, &u.q); }
    long p = 0, q = 0;                                    // u = p + q·sqrt(m) with small integers
    if (found) found = mpz_cmp_ui(mpq_denref(u.p.val), 1) == 0 && mpz_cmp_ui(mpq_denref(u.q.val), 1) == 0
                    && mpz_cmpabs_ui(mpq_numref(u.p.val), 2) <= 0 && mpz_cmpabs_ui(mpq_numref(u.q.val), 1) <= 0;
    if (found) { p = mpz_get_si(mpq_numref(u.p.val)); q = mpz_get_si(mpq_numref(u.q.val)); }
    long n = 0, d = 0;
    if (u.m == 2 && q == 1 && labs(p) == 1) { n = p < 0 ? 1 : 3; d = 8; }
    if (u.m == 3 && p == 2 && labs(q) == 1) { n = q < 0 ? 1 : 5; d = 12; }
    found = found && d;
    if (found) rat_set_si(turns, sign * n, d);
    qs_clear(&u);
    return found;
}

// atan(a), exact where tabulated, with atan(-u) = -atan(u); takes a
static Expr* fold_atan(Expr* a) {
    Rational t;
    rat_init(&t);
    bool known = atan_surd_value(a, &t);
    if (known) { free_expr(a); a = scaled(&t, create_pi()); }
    rat_clear(&t);
    if (known) return a;
    if (is_sum(a) && leading_negative(a))
        return fold_product(create_const(-1), fold_atan(fold_sum2(create_const(0), a, -1)));
    Rational c = rat_one_val(), turns;
    rat_init(&turns);
    Factors F = { NULL, 0 };
    gather(a, &c, &F);
    pair_sqrts(&c, &F);
    Expr* e;
    if (atan_value(&c, &F, &turns)) {
        free_factors(&F);
        e = scaled(&turns, create_pi());
    } else {
        bool neg = rat_is_negative(&c);
        if (neg) rat_neg(&c, &c);
        e = create_atan(make_product(&c, &F));
        if (neg) e = create_mul(create_const(-1), e);
    }
    free_expr(a);
    rat_clear(&c); rat_clear(&turns);
    return e;
}

// --- evaluation at a point: x0 = q, or x0 = q·π ---

// 1/r when r = c·sqrt(m_1)·...·sqrt(m_k), c ≠ 0: (1/(c·Π m_i))·Π sqrt(m_i); else NULL
static Expr* invert_surd(const Expr* r) {
    Rational c = rat_one_val();
    Factors F = { NULL, 0 };
    gather(r, &c, &F);
    bool ok = !rat_is_zero(&c);
    for (int i = 0; ok && i < F.n; i++) {
        ok = is_sqrt_const(F.f[i]) && !rat_is_zero(&F.f[i]->fn_arg->constant);
        if (ok) rat_mul(&c, &c, &F.f[i]->fn_arg->constant);
    }
    Expr* inv = NULL;
    if (ok) { mpq_inv(c.val, c.val); inv = make_product(&c, &F); }
    else free_factors(&F);
    rat_clear(&c);
    return inv;
}

// l / r folded when r is a rational times square roots, else kept; takes both
static Expr* fold_div(Expr* l, Expr* r) {
    Expr* inv = invert_surd(r);
    if (!inv) return create_div(l, r);
    free_expr(r);
    return fold_product(l, inv);
}

static Expr* fold_pow(Expr* l, Expr* r) {
    bool integer = is_const(r) && mpz_cmp_ui(mpq_denref(r->constant.val), 1) == 0;
    if (!is_const(l) || !integer) return create_pow(l, r);
    Rational v;
    rat_init(&v);
    rat_pow_int(&v, &l->constant, (int)mpz_get_si(mpq_numref(r->constant.val)));
    free_expr(l); free_expr(r);
    Expr* c = create_const_rat(&v);
    rat_clear(&v);
    return c;
}

static Expr* fold_binop(BinOp op, Expr* l, Expr* r) {
    switch (op) {
        case OP_ADD: return fold_sum2(l, r, 1);
        case OP_SUB: return fold_sum2(l, r, -1);
        case OP_MUL: return fold_product(l, r);
        case OP_DIV: return fold_div(l, r);
        case OP_POW: return fold_pow(l, r);
    }
    return NULL;
}


typedef struct {
    const Rational* q;
    bool pi;
} Point;

static Expr* rat_expr(const Rational* r) { return create_const_rat(r); }

// π^k (k >= 1): π itself lands in the sum's multiple-of-π slot, higher powers are factors
static Expr* pi_power(long k) { return k == 1 ? create_pi() : create_pow(create_pi(), create_const(k)); }

// p(x0) folded: Σ c_k·q^k·π^k at q·π (π is transcendental, so this is p's value in Q[π])
static Expr* poly_at(const Polynomial* p, const Point* x0) {
    Rational v;
    rat_init(&v);
    if (!x0->pi) { poly_eval_rat(&v, p, x0->q); Expr* c = rat_expr(&v); rat_clear(&v); return c; }
    Expr* acc = create_const(0);
    for (const Term* t = p->head; t; t = t->next) {
        rat_pow_int(&v, x0->q, t->exponent);
        rat_mul(&v, &v, &t->coeff);
        acc = fold_sum2(acc, t->exponent ? fold_product(rat_expr(&v), pi_power(t->exponent)) : rat_expr(&v), 1);
    }
    rat_clear(&v);
    return acc;
}

// p(x0) = a + r·π for p of degree <= 1 (r = 0 at a rational x0); false for higher degree at q·π
typedef struct {
    Rational a, r;
} Linear;

static bool linear_at(const Polynomial* p, const Point* x0, Linear* L) {
    if (!x0->pi) { poly_eval_rat(&L->a, p, x0->q); return true; }
    if (polynomial_degree(p) > 1) return false;
    poly_coeff(&L->a, p, 0);
    poly_coeff(&L->r, p, 1);
    rat_mul(&L->r, &L->r, x0->q);
    return true;
}

static Expr* at(const Expr* e, const Point* x0);

// exp(c): exp(0) = 1, else exp(c) itself
static Expr* fold_exp(Expr* c) {
    if (is_const(c) && rat_is_zero(&c->constant)) { free_expr(c); return create_const(1); }
    return create_exp(c);
}

// sin/cos(r) for rational r: sin(0) = 0, cos(0) = 1, sin(-r) = -sin(r), cos(-r) = cos(r)
static Expr* trig_of_rational(bool sine, const Rational* r) {
    if (rat_is_zero(r)) return create_const(sine ? 0 : 1);
    Rational a;
    rat_init(&a);
    rat_abs(&a, r);
    Expr* t = sine ? create_sin(rat_expr(&a)) : create_cos(rat_expr(&a));
    rat_clear(&a);
    return sine && rat_is_negative(r) ? fold_product(create_const(-1), t) : t;
}

// sin/cos(a + r·π) when 2r is an integer: a quarter-turn shift of sin/cos(a), sin(a + π/2) = cos(a)
static Expr* shifted_trig(bool sine, const Linear* L) {
    Rational t;
    rat_init(&t);
    rat_add(&t, &L->r, &L->r);
    mpz_t k;
    mpz_init(k);
    mpz_fdiv_r_ui(k, mpq_numref(t.val), 4);                // quarter turns mod 4
    int q = (int)mpz_get_ui(k);
    mpz_clear(k);
    rat_clear(&t);
    bool use_sin = sine ? q % 2 == 0 : q % 2 == 1;         // sin(a + π/2) = cos a, cos(a + π/2) = -sin a
    bool neg = sine ? q >= 2 : (q == 1 || q == 2);
    Expr* v = trig_of_rational(use_sin, &L->a);
    return neg ? fold_product(create_const(-1), v) : v;
}

// sin/cos of a polynomial argument at x0: exact at rational multiples of π, shifted by
// quarter turns around a rational, else symbolic
static Expr* fold_trig(const Expr* e, const Point* x0) {
    bool sine = e->type == EXPR_SIN;
    Linear L;
    rat_init(&L.a); rat_init(&L.r);
    Expr* r = NULL;
    if (e->fn_arg->type == EXPR_POLY && linear_at(e->fn_arg->poly, x0, &L)) {
        Rational t;
        rat_init(&t);
        rat_add(&t, &L.r, &L.r);
        bool quarter = mpz_cmp_ui(mpq_denref(t.val), 1) == 0;
        rat_clear(&t);
        if (rat_is_zero(&L.a) && !rat_is_zero(&L.r)) r = trig_value(sine ? TRIG_SIN : TRIG_COS, &L.r);
        else if (quarter) r = shifted_trig(sine, &L);
    }
    if (!r) {
        Expr* a = at(e->fn_arg, x0);
        if (is_const(a)) { r = trig_of_rational(sine, &a->constant); free_expr(a); }
        else r = sine ? create_sin(a) : create_cos(a);
    }
    rat_clear(&L.a); rat_clear(&L.r);
    return fold_product(create_const(1), r);
}

static Expr* at(const Expr* e, const Point* x0) {
    switch (e->type) {
        case EXPR_VAR: return x0->pi ? fold_product(rat_expr(x0->q), create_pi()) : rat_expr(x0->q);
        case EXPR_POLY: return poly_at(e->poly, x0);
        case EXPR_RATIONAL_FN: return fold_div(poly_at(e->ratfn.num, x0), poly_at(e->ratfn.den, x0));
        case EXPR_BINOP: return fold_binop(e->binop.op, at(e->binop.left, x0), at(e->binop.right, x0));
        case EXPR_NEG: return fold_product(create_const(-1), at(e->operand, x0));
        case EXPR_LN: return fold_ln(at(e->ln_arg, x0));
        case EXPR_ATAN: return fold_atan(at(e->fn_arg, x0));
        case EXPR_EXP: return fold_exp(at(e->fn_arg, x0));
        case EXPR_SIN: case EXPR_COS: return fold_trig(e, x0);
        default: return expr_copy(e);                     // CONST, SQRT, TRIG: no x inside
    }
}

Expr* expr_at(const Expr* e, const Rational* x0) { return at(e, &(Point){ x0, false }); }
Expr* expr_at_pi(const Expr* e, const Rational* q) { return at(e, &(Point){ q, true }); }

Expr* const_fold_sum(Expr** terms, int n) {
    Sum S = sum_new();
    for (int i = 0; i < n; i++) {
        add_all(&S, terms[i], 1);
        free_expr(terms[i]);
    }
    return sum_expr(&S, false);
}

// --- leading behaviour as x → ∞ ---

typedef struct {
    Expr* c; // e ~ c·x^deg, c a folded nonzero constant; NULL when not found
    int deg;
} Lead;

static Lead lead_of(const Expr* e);

static Lead lead_poly(const Polynomial* p) {
    if (poly_is_zero(p)) return (Lead){ NULL, 0 };
    Rational c;
    rat_init(&c);
    poly_coeff(&c, p, polynomial_degree(p));
    Lead L = { create_const_rat(&c), polynomial_degree(p) };
    rat_clear(&c);
    return L;
}

static bool is_zero_const(const Expr* e) { return is_const(e) && rat_is_zero(&e->constant); }

// the higher-degree side of l ± r; equal degrees add, failing when the leading terms cancel
static Lead lead_sum(Lead l, Lead r, int sign) {
    if (l.deg != r.deg) {
        Lead hi = l.deg > r.deg ? l : r, lo = l.deg > r.deg ? r : l;
        free_expr(lo.c);
        if (hi.c == r.c && sign < 0) hi.c = fold_product(create_const(-1), hi.c);
        return hi;
    }
    Lead s = { fold_sum2(l.c, r.c, sign), l.deg };
    if (is_zero_const(s.c)) { free_expr(s.c); s.c = NULL; }
    return s;
}

// x-power base^n for an integer constant n
static Lead lead_pow(Lead b, const Expr* n) {
    bool integer = is_const(n) && mpz_cmp_ui(mpq_denref(n->constant.val), 1) == 0;
    if (!integer) { free_expr(b.c); return (Lead){ NULL, 0 }; }
    long k = mpz_get_si(mpq_numref(n->constant.val));
    return (Lead){ fold_pow(b.c, expr_copy(n)), b.deg * (int)k };
}

static Lead lead_binop(const Expr* e) {
    Lead l = lead_of(e->binop.left);
    if (!l.c) return l;
    if (e->binop.op == OP_POW) return lead_pow(l, e->binop.right);
    Lead r = lead_of(e->binop.right);
    if (!r.c) { free_expr(l.c); return r; }
    switch (e->binop.op) {
        case OP_ADD: return lead_sum(l, r, 1);
        case OP_SUB: return lead_sum(l, r, -1);
        case OP_MUL: return (Lead){ fold_product(l.c, r.c), l.deg + r.deg };
        default:     return (Lead){ fold_div(l.c, r.c), l.deg - r.deg };
    }
}

static Lead lead_ratfn(const Expr* e) {
    Lead n = lead_poly(e->ratfn.num), d = lead_poly(e->ratfn.den);
    if (n.c && d.c) return (Lead){ fold_div(n.c, d.c), n.deg - d.deg };
    free_expr(n.c); free_expr(d.c);
    return (Lead){ NULL, 0 };
}

// e contains x
static bool has_x(const Expr* e) {
    switch (e->type) {
        case EXPR_VAR: return true;
        case EXPR_POLY: return polynomial_degree(e->poly) > 0;
        case EXPR_RATIONAL_FN: return true;
        case EXPR_BINOP: return has_x(e->binop.left) || has_x(e->binop.right);
        case EXPR_NEG: return has_x(e->operand);
        case EXPR_LN: return has_x(e->ln_arg);
        case EXPR_ATAN: case EXPR_SQRT: case EXPR_EXP: case EXPR_SIN: case EXPR_COS: return has_x(e->fn_arg);
        default: return e->type == EXPR_ROOTSUM;
    }
}

static Lead lead_of(const Expr* e) {
    if (!has_x(e)) {                                      // a constant: itself, folded
        Rational zero = rat_zero_val();
        Lead L = { expr_at(e, &zero), 0 };
        rat_clear(&zero);
        if (is_zero_const(L.c)) { free_expr(L.c); L.c = NULL; }
        return L;
    }
    switch (e->type) {
        case EXPR_VAR: return (Lead){ create_const(1), 1 };
        case EXPR_POLY: return lead_poly(e->poly);
        case EXPR_RATIONAL_FN: return lead_ratfn(e);
        case EXPR_BINOP: return lead_binop(e);
        case EXPR_NEG: {
            Lead L = lead_of(e->operand);
            if (L.c) L.c = fold_product(create_const(-1), L.c);
            return L;
        }
        default: return (Lead){ NULL, 0 };              // ln, atan, RootSum of x: not polynomial-like
    }
}

Expr* expr_lead(const Expr* e, int* deg) {
    Lead L = lead_of(e);
    *deg = L.deg;
    return L.c;
}
