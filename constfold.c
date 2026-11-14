// Exact folding of constant expressions produced by evaluating an antiderivative at a
// rational point. A term is c · f_1 · ... · f_k with c rational and f_i non-rational
// factors (radicals, sqrt, trig constants, ln, atan, π); a sum is a rational part, a
// multiple of π, logs of rationals c·ln(v) merged by v, and other terms merged when their
// factors print the same.
#include "constfold.h"
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

// sqrt(m)·sqrt(m) = m
static void pair_sqrts(Rational* c, Factors* F) {
    for (int i = 0; i < F->n; i++)
        for (int j = i + 1; F->f[i] && is_sqrt_const(F->f[i]) && j < F->n; j++) {
            if (!F->f[j] || !is_sqrt_const(F->f[j]) || !rat_eq(&F->f[i]->fn_arg->constant, &F->f[j]->fn_arg->constant)) continue;
            rat_mul(c, c, &F->f[i]->fn_arg->constant);
            free_expr(F->f[i]); free_expr(F->f[j]);
            F->f[i] = F->f[j] = NULL;
        }
    int k = 0;
    for (int i = 0; i < F->n; i++) if (F->f[i]) F->f[k++] = F->f[i];
    F->n = k;
}

// c · Π F as an expression; takes F's factors
static Expr* make_product(const Rational* c, Factors* F) {
    if (rat_is_zero(c)) { free_factors(F); return create_const(0); }
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

// the folded sum; the rational part goes first when rational_first, else last
static Expr* sum_expr(Sum* S, bool rational_first) {
    merge_logs(S);
    Expr* e = NULL;
    bool rat = !rat_is_zero(&S->rational);
    if (rat && rational_first) e = create_const_rat(&S->rational);
    for (int i = 0; i < S->nt; i++) {
        if (!rat_is_zero(&S->t[i].c)) e = append(e, scaled(&S->t[i].c, S->t[i].body));
        else free_expr(S->t[i].body);
        free(S->t[i].key);
        rat_clear(&S->t[i].c);
    }
    for (int i = 0; i < S->nl; i++) {
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

// atan(a), exact where tabulated, with atan(-u) = -atan(u); takes a
static Expr* fold_atan(Expr* a) {
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

// --- evaluation at a point ---

static Expr* const_value(const Polynomial* p, const Rational* x0) {
    Rational v;
    rat_init(&v);
    poly_eval_rat(&v, p, x0);
    Expr* c = create_const_rat(&v);
    rat_clear(&v);
    return c;
}

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

Expr* expr_at(const Expr* e, const Rational* x0) {
    switch (e->type) {
        case EXPR_VAR: return create_const_rat(x0);
        case EXPR_POLY: return const_value(e->poly, x0);
        case EXPR_RATIONAL_FN: return fold_div(const_value(e->ratfn.num, x0), const_value(e->ratfn.den, x0));
        case EXPR_BINOP: return fold_binop(e->binop.op, expr_at(e->binop.left, x0), expr_at(e->binop.right, x0));
        case EXPR_NEG: return fold_product(create_const(-1), expr_at(e->operand, x0));
        case EXPR_LN: return fold_ln(expr_at(e->ln_arg, x0));
        case EXPR_ATAN: return fold_atan(expr_at(e->fn_arg, x0));
        default: return expr_copy(e);                     // CONST, SQRT, TRIG: no x inside
    }
}

Expr* const_fold_sum(Expr** terms, int n) {
    Sum S = sum_new();
    for (int i = 0; i < n; i++) {
        add_all(&S, terms[i], 1);
        free_expr(terms[i]);
    }
    return sum_expr(&S, false);
}
