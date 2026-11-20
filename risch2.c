// Two-level Risch integration: θ2 = ln(u) over K1 = Q(x)(θ1) (risch2.h).
//
// For f = Σ_{k=0}^{n} a_k·θ2^k with a_k ∈ K1, Liouville's theorem gives (as for one ln over Q(x),
// risch.c) ∫ f = Σ_{k=1}^{n+1} b_k·θ2^k + ∫ (level 0), b_{n+1} constant, b_k ∈ K1, and
//   a_k = D(b_k) + (k+1)·b_{k+1}·η,  η = D(u)/u.
// From the top, g_k = a_k - (k+1)·r_{k+1}·η must satisfy g_k = D(r_k) + c·η with r_k ∈ K1 and a
// constant c (limited integration over K1); then b_{k+1} = r_{k+1} + c/(k+1). Level 0 is any
// integral over K1 (the one-tower integrator).
//
// Limited integration over K1 uses only Hermite reduction: g = D(R) + h/d + P (d square-free,
// P a polynomial in θ1, Laurent for exp). D(R') for R' ∈ K1 has no simple normal poles, so the
// simple part h/d of an element of D(K1) + Q·η is c times that of η; P_η lies in Q(x) (D(u)/u
// has no θ1^k part, k ≠ 0, for exp, and is proper in θ1 for ln). So g = D(r) + c·η iff
//   1. h_g/d_g = c·h_η/d_η;
//   2. the θ1^k parts of P_g, k ≠ 0, integrate in K1 (exp: an RDE per k; ln: the levels of
//      risch.c, which leave θ1^0 part q_g = p_0 - b_1·L1);
//   3. at θ1^0, q_g - c·p_η has an integral in Q(x), plus d·θ1 for ln (d·L1 = D(d·θ1)): the
//      Hermite simple parts satisfy s(q_g) = c·s(p_η) + d·L1.
// Every answer is checked exactly: g = D(r) + c·η, and D(Σ b_k θ2^k) + g_0 = f. θ2 must be
// transcendental over K1: η ∉ D(K1), i.e. not (h_η = 0 and s(p_η) ∈ Q·L1).
#include "risch2.h"
#include "k2poly.h"
#include "factor.h"
#include "rischk1.h"
#include "integrate.h"
#include "matrix.h"
#include "expsum.h"
#include <stdlib.h>
#include <string.h>

const char TOWER2_NOT_APPLICABLE[] = "ln over another exp or ln tower is supported only for polynomials in the outer ln for now";

// --- expression -> K1[θ2] ---

typedef struct {
    const Tower* t;
    const char* top;    // the printed argument of the outer ln
    const char* why;
} Conv2;

static bool is_top(const Expr* e, const char* top) {
    if (e->type != EXPR_LN) return false;
    char* s = expr_to_string(e->ln_arg);
    bool same = strcmp(s, top) == 0;
    free(s);
    return same;
}

static bool has_top(const Expr* e, const char* top) {
    if (is_top(e, top)) return true;
    switch (e->type) {
        case EXPR_BINOP: return has_top(e->binop.left, top) || has_top(e->binop.right, top);
        case EXPR_NEG: return has_top(e->operand, top);
        case EXPR_LN: return has_top(e->ln_arg, top);
        case EXPR_EXP: case EXPR_SIN: case EXPR_COS: case EXPR_ATAN: case EXPR_SQRT: return has_top(e->fn_arg, top);
        default: return false;
    }
}

static K2Frac conv2(const Expr* e, Conv2* C);

// e ∈ K1 (no outer ln)
static bool conv_k1(const Expr* e, Conv2* C, KFrac* out) {
    C->why = exp_form_in(e, C->t, out);
    return !C->why;
}

static K2Poly k2_one(TowerKind kind) { return k2_const(kfrac_make(kp_one(), kp_one(), kind), kind); }

// p/1 (takes p)
static K2Frac k2frac_poly(K2Poly p) { return (K2Frac){ p, k2_one(p.kind) }; }

static K2Frac k2frac_zero(TowerKind kind) { return k2frac_poly(k2_zeros(0, kind)); }

// b^k for an integer k (k < 0: the inverse first)
static K2Frac conv2_pow(const Expr* e, Conv2* C) {
    TowerKind kind = C->t->kind;
    long k = 0;
    if (!small_integer(e->binop.right, &k)) { C->why = "exponents must be small integers"; return k2frac_zero(kind); }
    K2Frac b = conv2(e->binop.left, C), r = k2frac_poly(k2_one(kind));
    if (!C->why && k < 0) {
        if (b.n.n == 0) C->why = "division by zero";
        else { K2Frac inv = k2frac_inv(&b); k2frac_free(&b); b = inv; k = -k; }
    }
    for (long i = 0; i < k && !C->why; i++) { K2Frac m = k2frac_mul(&r, &b); k2frac_free(&r); r = m; }
    k2frac_free(&b);
    return r;
}

static K2Frac conv2_binop(const Expr* e, Conv2* C) {
    if (e->binop.op == OP_POW) return conv2_pow(e, C);
    K2Frac l = conv2(e->binop.left, C), r = conv2(e->binop.right, C), out;
    if (C->why) out = k2frac_zero(C->t->kind);
    else if (e->binop.op == OP_MUL) out = k2frac_mul(&l, &r);
    else if (e->binop.op == OP_DIV) {
        if (r.n.n == 0) { C->why = "division by zero"; out = k2frac_zero(C->t->kind); }
        else if (k2_deg(&r.n) == 0 && k2_deg(&r.d) == 0) {  // by an element of K1: scale the numerator
            KFrac ic = kfrac_inv(&r.n.a[0], C->t->kind), dc = kfrac_mul(&ic, &r.d.a[0], C->t->kind);
            K2Poly s = k2_const(dc, C->t->kind);
            out = (K2Frac){ k2_mul(&l.n, &s), k2_copy(&l.d) };
            k2_free(&s); kfrac_free(&ic);
        } else { K2Frac ir = k2frac_inv(&r); out = k2frac_mul(&l, &ir); k2frac_free(&ir); }
    } else out = k2frac_add(&l, &r, e->binop.op == OP_ADD ? 1 : -1);
    k2frac_free(&l); k2frac_free(&r);
    return out;
}

static K2Frac conv2(const Expr* e, Conv2* C) {
    TowerKind kind = C->t->kind;
    if (C->why) return k2frac_zero(kind);
    if (is_top(e, C->top)) {                              // θ2
        K2Poly p = k2_zeros(0, kind);
        k2_accumulate(&p, 1, kfrac_make(kp_one(), kp_one(), kind));
        return k2frac_poly(p);
    }
    if (!has_top(e, C->top)) {
        KFrac k;
        return conv_k1(e, C, &k) ? k2frac_poly(k2_const(k, kind)) : k2frac_zero(kind);
    }
    switch (e->type) {
        case EXPR_BINOP: return conv2_binop(e, C);
        case EXPR_NEG: {
            K2Frac z = k2frac_zero(kind), f = conv2(e->operand, C), r = k2frac_add(&z, &f, -1);
            k2frac_free(&z); k2frac_free(&f);
            return r;
        }
        default: C->why = "the outer ln inside another function is not supported yet"; return k2frac_zero(kind);
    }
}

// e with every outer ln(top) replaced by x: what remains decides the tower K1
static Expr* without_top(const Expr* e, const char* top) {
    if (is_top(e, top)) return create_var("x");
    switch (e->type) {
        case EXPR_BINOP: {
            Expr *l = without_top(e->binop.left, top), *r = without_top(e->binop.right, top);
            switch (e->binop.op) {
                case OP_ADD: return create_add(l, r);
                case OP_SUB: return create_sub(l, r);
                case OP_MUL: return create_mul(l, r);
                case OP_DIV: return create_div(l, r);
                default: return create_pow(l, r);
            }
        }
        case EXPR_NEG: return create_neg(without_top(e->operand, top));
        default: return expr_copy(e);
    }
}

// the distinct ln arguments of e, nested ones (containing exp or ln) first
typedef struct {
    const Expr** arg;
    int n;
} LnArgs;

static void collect_ln(const Expr* e, LnArgs* A) {
    switch (e->type) {
        case EXPR_LN: {
            char* s = expr_to_string(e->ln_arg);
            bool seen = false;
            for (int i = 0; i < A->n && !seen; i++) { char* t = expr_to_string(A->arg[i]); seen = strcmp(s, t) == 0; free(t); }
            free(s);
            if (!seen) { A->arg = realloc(A->arg, (size_t)(A->n + 1) * sizeof(Expr*)); A->arg[A->n++] = e->ln_arg; }
            collect_ln(e->ln_arg, A);
            return;
        }
        case EXPR_BINOP: collect_ln(e->binop.left, A); collect_ln(e->binop.right, A); return;
        case EXPR_NEG: collect_ln(e->operand, A); return;
        case EXPR_EXP: case EXPR_SIN: case EXPR_COS: case EXPR_ATAN: case EXPR_SQRT: collect_ln(e->fn_arg, A); return;
        default: return;
    }
}

static bool nested(const Expr* a) { return expr_has_exp(a) || expr_has_ln(a); }

// f in two levels with outer ln(u)
typedef struct {
    Tower t1;
    KFrac u;
    const Expr* u_expr;   // points into f
    K2Frac f;
    const char* why;      // when there is no such form: the reason from a tower that fit
} Form2;

// k is a rational constant
static bool kf_is_constant(const KFrac* k) {
    Rational c;
    rat_init(&c);
    bool yes = kf_is_zero(k) || (kp_is_one(&k->d) && k->n.n == 1 && k->n.lo == 0 && qfn_constant(&k->n.a[0], &c));
    rat_clear(&c);
    return yes;
}

// f with outer ln(A), or false (out->why set when the towers fit but f's shape does not): K1
// is the tower of the rest of f together with A
static bool form_with_top(const Expr* f, const Expr* A, Form2* out) {
    char* top = expr_to_string(A);
    Expr* rest = create_add(without_top(f, top), expr_copy(A));
    ExpForm F = exp_form(rest);
    free_expr(rest);
    bool ok = !F.unsupported;
    if (ok) {
        out->t1 = tower_copy(&F.t);
        out->u_expr = A;
        Conv2 C = { &out->t1, top, NULL };
        ok = conv_k1(A, &C, &out->u);
        if (ok && kf_is_constant(&out->u)) { kfrac_free(&out->u); ok = false; }   // ln of a constant
        if (ok) {
            out->f = conv2(f, &C);
            if (C.why) { k2frac_free(&out->f); kfrac_free(&out->u); ok = false; out->why = C.why; }
        }
        if (!ok) tower_free(&out->t1);
    }
    exp_form_free(&F);
    free(top);
    return ok;
}

// false with out->why: the reason from a tower that fit, else TOWER2_NOT_APPLICABLE
static bool form2_of(const Expr* f, Form2* out) {
    LnArgs A = { NULL, 0 };
    collect_ln(f, &A);
    bool ok = false;
    out->why = TOWER2_NOT_APPLICABLE;
    for (int pass = 0; pass < 2 && !ok; pass++)            // nested arguments first, then the rest from the last
        for (int i = A.n - 1; i >= 0 && !ok; i--)
            if (nested(A.arg[i]) == (pass == 0)) ok = form_with_top(f, A.arg[i], out);
    free(A.arg);
    return ok;
}

static void form2_free(Form2* F) {
    k2frac_free(&F->f);
    kfrac_free(&F->u);
    tower_free(&F->t1);
}

// --- limited integration over K1 ---

// g = D(R) + S + P: R and S (= h/d, d square-free) from Hermite reduction, P a polynomial in θ1
typedef struct {
    KFrac R, S;
    LPoly P;
} Split;

static void split_free(Split* s) { kfrac_free(&s->R); kfrac_free(&s->S); lpoly_free(&s->P); }

static bool split_k1(const KFrac* g, const Tower* t, Split* out) {
    TowerKind kind = t->kind;
    if (kp_is_one(&g->d)) { *out = (Split){ kf_zero(), kf_zero(), kp_copy(&g->n) }; return true; }
    ExpForm F = { *t, kp_copy(&g->n), kp_copy(&g->d), NULL };  // shares t
    LPoly b = proper_part(&F);
    Reduced Rd = hermite_reduce(&b, &F);
    out->R = kf_copy(&Rd.g);
    out->S = Rd.h.n ? kfrac_make(kp_copy(&Rd.h), kp_copy(&Rd.d), kind) : kf_zero();
    KFrac DR = kfrac_derivative(&out->R, t), a = kfrac_sub(g, &DR, kind), rem = kfrac_sub(&a, &out->S, kind);
    bool ok = kp_is_one(&rem.d);
    out->P = ok ? kp_copy(&rem.n) : (LPoly){ 0, 0, NULL };
    kfrac_free(&DR); kfrac_free(&a); kfrac_free(&rem);
    reduced_free(&Rd); lpoly_free(&b); lpoly_free(&F.f); lpoly_free(&F.den);
    if (!ok) { kfrac_free(&out->R); kfrac_free(&out->S); }
    return ok;
}

// a / b
static KFrac kf_div(const KFrac* a, const KFrac* b, TowerKind kind) {
    KFrac ib = kfrac_inv(b, kind), q = kfrac_mul(a, &ib, kind);
    kfrac_free(&ib);
    return q;
}

// q is a rational constant c
static bool kf_constant(const KFrac* q, Rational* c) {
    return kp_is_one(&q->d) && q->n.n == 1 && q->n.lo == 0 && qfn_constant(&q->n.a[0], c);
}

// Σ x_i·basis_i = target in Q(x) for constants x_i
typedef struct {
    const QFn* target;
    const QFn* basis;
    int m;
} ConstSystem;

// UNDETERMINED (dependent basis) never ends in a proof: the caller reports it as unsupported
typedef enum { SOLVED, NO_SOLUTION, UNDETERMINED } Solve;

// the numerators of the system over one common denominator, coefficient by coefficient
static Solve solve_constants(const ConstSystem* S, Rational* x) {
    int m = S->m;
    Polynomial* D = poly_copy(S->target->den);
    for (int i = 0; i < m; i++) { Polynomial* p = poly_multiply(D, S->basis[i].den); free_polynomial(D); D = p; }
    Polynomial** N = malloc((size_t)(m + 1) * sizeof(Polynomial*));
    int deg = 0;
    for (int i = 0; i <= m; i++) {
        const QFn* q = i < m ? &S->basis[i] : S->target;
        Polynomial* k = poly_divide(D, q->den);
        N[i] = poly_multiply(q->num, k);
        free_polynomial(k);
        if (!poly_is_zero(N[i]) && polynomial_degree(N[i]) > deg) deg = polynomial_degree(N[i]);
    }
    Solve st = poly_is_zero(N[m]) ? SOLVED : NO_SOLUTION;
    if (m > 0) {
        Matrix* A = matrix_create(deg + 1, m);
        Vector* b = vector_create(deg + 1);
        for (int r = 0; r <= deg; r++) {
            for (int i = 0; i < m; i++) poly_coeff(&MAT_AT(A, r, i), N[i], r);
            poly_coeff(&b->data[r], N[m], r);
        }
        LinSysResult* L = linsys_solve(A, b);
        st = L->status == LINSYS_INCONSISTENT ? NO_SOLUTION : L->status == LINSYS_INFINITE ? UNDETERMINED : SOLVED;
        for (int i = 0; st == SOLVED && i < m; i++) rat_set(&x[i], &L->solution->data[i]);
        linsys_result_free(L);
        matrix_free(A); vector_free(b);
    }
    for (int i = 0; i <= m; i++) free_polynomial(N[i]);
    free(N);
    free_polynomial(D);
    return st;
}

// Q(x) as a constant-coefficient element of K1
static KFrac kf_of_qfn(const QFn* q, TowerKind kind) {
    return kfrac_make(lpoly_term(qfn_copy(q), 0), kp_one(), kind);
}

// η and the tower, split once: η's Hermite parts, its θ1^0 polynomial part p (in Q(x)) and p's
// Hermite split
typedef struct {
    const Tower* t;
    KFrac eta;
    Split sp;
    QFn p, s, R0;     // p = s + D(R0), s simple
    QFn L1;           // ln base: L1 = D(θ1); else 0
} Eta;

typedef enum { LIM_OK, LIM_NONE, LIM_UNSUPPORTED } LimStatus;

// q = s + D(R0) over Q(x), s the Hermite simple part; false (s = R0 = 0) when Hermite's check fails
static bool hermite_q(const QFn* q, QFn* s, QFn* R0) {
    HermiteSplit h;
    bool ok = hermite_split(q->num, q->den, &h);
    *s = ok ? qfn_from(poly_copy(h.a), poly_copy(h.dstar)) : qfn_zero();
    *R0 = ok ? qfn_from(poly_copy(h.rnum), poly_copy(h.rden)) : qfn_zero();
    hermite_split_free(&h);
    return ok;
}

static bool eta_init(const KFrac* u, const Tower* t, Eta* E) {
    E->t = t;
    KFrac Du = kfrac_derivative(u, t);
    E->eta = kf_div(&Du, u, t->kind);
    kfrac_free(&Du);
    E->L1 = t->kind == TOWER_LOG ? tower_L(t) : qfn_zero();
    bool ok = split_k1(&E->eta, t, &E->sp);
    if (ok) ok = E->sp.P.n == 0 || (E->sp.P.n == 1 && E->sp.P.lo == 0);  // P_η ∈ Q(x)
    E->p = ok ? kp_coeff(&E->sp.P, 0) : qfn_zero();
    return hermite_q(&E->p, &E->s, &E->R0) && ok;
}

static void eta_free(Eta* E) {
    kfrac_free(&E->eta); split_free(&E->sp);
    qfn_free(&E->p); qfn_free(&E->s); qfn_free(&E->R0); qfn_free(&E->L1);
}

// θ2 = ln(u) may lie in K1 (η ∈ D(K1)): h_η = 0 and s_η ∈ Q·L1 (undetermined counts as "may")
static bool eta_in_dk1(const Eta* E) {
    if (!kf_is_zero(&E->sp.S)) return false;
    if (E->t->kind != TOWER_LOG) return qfn_is_zero(&E->s);
    Rational d;
    rat_init(&d);
    Solve st = solve_constants(&(ConstSystem){ &E->s, &E->L1, 1 }, &d);
    rat_clear(&d);
    return st != NO_SOLUTION;
}

// the θ1^k parts (k ≠ 0) of a polynomial P integrated in K1: D(Y) equals them, and q is the θ1^0
// part left over (for ln, q = p_0 - b_1·L1 from the levels)
typedef struct {
    LPoly Y;
    QFn q;
} PolyPart;

static LimStatus poly_part(const LPoly* P, const Eta* E, PolyPart* out) {
    const Tower* t = E->t;
    ExpForm F = { *t, kp_copy(P), kp_one(), NULL };       // shares t
    LimStatus st = LIM_OK;
    out->Y = (LPoly){ 0, 0, NULL };
    out->q = qfn_zero();
    if (t->kind == TOWER_EXP) {
        int bad = 0;
        if (exp_part(&F, &out->Y, &bad)) { qfn_free(&out->q); out->q = kp_coeff(P, 0); }
        else { st = LIM_NONE; out->Y = (LPoly){ 0, 0, NULL }; }  // exp_part freed Y
    } else {
        LogSolution S = { NULL, qfn_zero(), P->n > 0 ? P->lo + P->n - 1 : 0 };
        S.b = malloc((size_t)(S.n + 2) * sizeof(QFn));
        for (int k = 0; k <= S.n + 1; k++) S.b[k] = qfn_zero();
        if (solve_log_levels(&F, &S)) { out->Y = log_part_poly(&S); qfn_free(&out->q); out->q = qfn_copy(&S.g0); }
        else st = LIM_NONE;
        for (int k = 0; k <= S.n + 1; k++) qfn_free(&S.b[k]);
        free(S.b);
        qfn_free(&S.g0);
    }
    lpoly_free(&F.f); lpoly_free(&F.den);
    return st;
}

// c (and d for ln) from s_g = c·s_η + d·L1, c fixed when known
typedef struct {
    Rational c, d;
    bool c_known;
} Consts;

static Solve solve_cd(const QFn* s_g, const Eta* E, Consts* K) {
    bool ln = E->t->kind == TOWER_LOG;
    QFn basis[2];
    int m = 0;
    QFn target = qfn_copy(s_g);
    if (K->c_known) { QFn cs = { poly_scale(E->s.num, &K->c), poly_copy(E->s.den) }, tt = qfn_add(&target, &cs, -1); qfn_free(&target); qfn_free(&cs); target = tt; }
    else basis[m++] = E->s;
    if (ln) basis[m++] = E->L1;
    Rational x[2];
    rat_init(&x[0]); rat_init(&x[1]);
    Solve st = solve_constants(&(ConstSystem){ &target, basis, m }, x);
    if (st == SOLVED && !K->c_known) { rat_set(&K->c, &x[0]); K->c_known = true; }
    if (st == SOLVED && ln) rat_set(&K->d, &x[m - 1]);
    rat_clear(&x[0]); rat_clear(&x[1]);
    qfn_free(&target);
    return st;
}

// g = D(r) + c·η with r ∈ K1, c ∈ Q
typedef struct {
    KFrac r;
    Rational c;   // initialized by the caller
} LimResult;

// what steps 1-3 find for g: the constants and g's polynomial part integrated
typedef struct {
    Consts K;
    PolyPart W;
} LimWork;

// steps 1-3 for g's split; LIM_NONE when no constants exist
static LimStatus lim_constants(const Split* sg, const Eta* E, LimWork* L) {
    Consts* K = &L->K;
    PolyPart* W = &L->W;
    LimStatus st = LIM_OK;
    if (kf_is_zero(&E->sp.S)) st = kf_is_zero(&sg->S) ? LIM_OK : LIM_NONE;   // 1. simple parts
    else if (kf_is_zero(&sg->S)) K->c_known = true;        // c = 0
    else {
        KFrac q = kf_div(&sg->S, &E->sp.S, E->t->kind);
        st = kf_constant(&q, &K->c) ? LIM_OK : LIM_NONE;
        K->c_known = true;
        kfrac_free(&q);
    }
    *W = (PolyPart){ { 0, 0, NULL }, qfn_zero() };
    if (st == LIM_OK) { qfn_free(&W->q); st = poly_part(&sg->P, E, W); }             // 2.
    QFn s, R0;
    if (st == LIM_OK) {                                                               // 3.
        Solve sv = hermite_q(&W->q, &s, &R0) ? solve_cd(&s, E, K) : UNDETERMINED;
        st = sv == SOLVED ? LIM_OK : sv == NO_SOLUTION ? LIM_NONE : LIM_UNSUPPORTED;
        qfn_free(&s); qfn_free(&R0);
    }
    return st;
}

// g, its split and what steps 1-3 found
typedef struct {
    const KFrac* g;
    const Split* sg;
    const LimWork* L;
} LimParts;

// r = R_g - c·R_η + Y + d·θ1 + ∫(q - c·p_η - d·L1), the last rational by step 3; false when
// that or the exact check g = D(r) + c·η fails
static bool lim_answer(const LimParts* X, const Eta* E, KFrac* r) {
    TowerKind kind = E->t->kind;
    const KFrac* g = X->g;
    const Split* sg = X->sg;
    const Consts* K = &X->L->K;
    const PolyPart* W = &X->L->W;
    QFn cp = { poly_scale(E->p.num, &K->c), poly_copy(E->p.den) }, dl = { poly_scale(E->L1.num, &K->d), poly_copy(E->L1.den) };
    QFn q1 = qfn_add(&W->q, &cp, -1), q2 = qfn_add(&q1, &dl, -1), s2, R2;
    bool ok = hermite_q(&q2, &s2, &R2) && qfn_is_zero(&s2);
    KFrac cR = kf_scaled(&E->sp.R, &K->c, kind), a = kfrac_sub(&sg->R, &cR, kind);
    KFrac y = kfrac_make(kp_copy(&W->Y), kp_one(), kind), b = kfrac_add(&a, &y, kind);
    KFrac dth = kfrac_make(lpoly_term(qfn_rat(&K->d), 1), kp_one(), kind), c1 = kfrac_add(&b, &dth, kind);
    KFrac r0 = kf_of_qfn(&R2, kind);
    *r = kfrac_add(&c1, &r0, kind);
    KFrac Dr = kfrac_derivative(r, E->t), ce = kf_scaled(&E->eta, &K->c, kind), e1 = kfrac_sub(g, &Dr, kind), e2 = kfrac_sub(&e1, &ce, kind);
    ok = ok && kf_is_zero(&e2);
    if (!ok) kfrac_free(r);
    qfn_free(&cp); qfn_free(&dl); qfn_free(&q1); qfn_free(&q2); qfn_free(&s2); qfn_free(&R2);
    kfrac_free(&cR); kfrac_free(&a); kfrac_free(&y); kfrac_free(&b); kfrac_free(&dth); kfrac_free(&c1); kfrac_free(&r0);
    kfrac_free(&Dr); kfrac_free(&ce); kfrac_free(&e1); kfrac_free(&e2);
    return ok;
}

static LimStatus lim_integrate(const KFrac* g, const Eta* E, LimResult* out) {
    Split sg;
    if (!split_k1(g, E->t, &sg)) return LIM_UNSUPPORTED;
    LimWork L;
    rat_init(&L.K.c); rat_init(&L.K.d);
    L.K.c_known = false;
    LimStatus st = lim_constants(&sg, E, &L);
    if (st == LIM_OK && !lim_answer(&(LimParts){ g, &sg, &L }, E, &out->r)) st = LIM_UNSUPPORTED;
    if (st == LIM_OK) rat_set(&out->c, &L.K.c);
    lpoly_free(&L.W.Y); qfn_free(&L.W.q);
    split_free(&sg);
    rat_clear(&L.K.c); rat_clear(&L.K.d);
    return st;
}

// --- the levels in θ2 ---

#define LEVEL_FAILS "has no elementary antiderivative (Risch: a level of the outer ln needs a logarithm other than it)"

// D(Σ b_k θ2^k) = Σ (D(b_k)·θ2^k + k·b_k·η·θ2^(k-1))
static K2Poly k2_derivative(const K2Poly* B, const Eta* E) {
    K2Poly r = k2_zeros(0, B->kind);
    for (int k = 0; k < B->n; k++) {
        k2_accumulate(&r, k, kfrac_derivative(&B->a[k], E->t));
        if (k == 0) continue;
        Rational kk = rat_from_int(k);
        KFrac be = kfrac_mul(&B->a[k], &E->eta, B->kind);
        k2_accumulate(&r, k - 1, kf_scaled(&be, &kk, B->kind));
        kfrac_free(&be);
        rat_clear(&kk);
    }
    k2_trim(&r);
    return r;
}

// ∫ f = Σ_{k>=1} b_k θ2^k + ∫ g0
typedef struct {
    K2Poly B;
    KFrac g0;
} Levels;

// level k's result: r_{k+1} from the level above and the c found at level k
typedef struct {
    int k;
    const KFrac* r_next;
    const Rational* c;
} Step;

// b_{k+1} = r_{k+1} + c/(k+1), with g_k = D(r_k) + c·η
static void set_level(Levels* L, const Step* st) {
    TowerKind kind = L->B.kind;
    int k = st->k;
    const KFrac* r_next = st->r_next;
    const Rational* c = st->c;
    Rational kk = rat_from_int(k + 1), d;
    rat_init(&d);
    rat_div(&d, c, &kk);
    QFn dq = qfn_rat(&d);
    KFrac dk = kf_of_qfn(&dq, kind);
    kfrac_free(&L->B.a[k + 1]);
    L->B.a[k + 1] = kfrac_add(r_next, &dk, kind);
    kfrac_free(&dk); qfn_free(&dq);
    rat_clear(&kk); rat_clear(&d);
}

// a polynomial P in θ2 to integrate, in the towers of F
typedef struct {
    const Form2* F;
    const Eta* E;
    const K2Poly* P;
} PolyIn;

// the levels from the top; LIM_NONE at one that fails
static LimStatus solve_levels(const PolyIn* in, Levels* L) {
    const Form2* F = in->F;
    const Eta* E = in->E;
    const K2Poly* P = in->P;
    TowerKind kind = F->t1.kind;
    int n = P->n - 1;                                     // -1 for P = 0
    L->B = k2_zeros(n + 2 < 2 ? 2 : n + 2, kind);
    L->g0 = kf_zero();
    KFrac r_next = kf_zero();
    LimStatus st = LIM_OK;
    for (int k = n; st == LIM_OK && k >= 1; k--) {
        Rational kk = rat_from_int(k + 1);
        KFrac re = kfrac_mul(&r_next, &E->eta, kind), t = kf_scaled(&re, &kk, kind), g = kfrac_sub(&P->a[k], &t, kind);
        LimResult R;
        rat_init(&R.c);
        st = lim_integrate(&g, E, &R);
        if (st == LIM_OK) { set_level(L, &(Step){ k, &r_next, &R.c }); kfrac_free(&r_next); r_next = R.r; }
        kfrac_free(&re); kfrac_free(&t); kfrac_free(&g);
        rat_clear(&kk); rat_clear(&R.c);
    }
    if (st == LIM_OK) {
        kfrac_free(&L->B.a[1]);
        L->B.a[1] = kf_copy(&r_next);
        KFrac re = kfrac_mul(&r_next, &E->eta, kind), zero = kf_zero();
        kfrac_free(&L->g0);
        L->g0 = kfrac_sub(n >= 0 ? &P->a[0] : &zero, &re, kind);
        kfrac_free(&re); kfrac_free(&zero);
    }
    kfrac_free(&r_next);
    k2_trim(&L->B);
    return st;
}

// b·ln(u)^k, printed as num·ln(u)^k / den when b is a fraction in θ1
static Expr* level_term(const KFrac* b, const Form2* F, int k) {
    Expr* ln = create_ln(expr_copy(F->u_expr));
    if (k > 1) ln = create_pow(ln, create_const(k));
    Rational c;
    rat_init(&c);
    bool constant = kf_constant(b, &c);
    Expr* e = constant ? create_mul(create_const_rat(&c), ln) : NULL;  // 1·ln and -1·ln print bare
    rat_clear(&c);
    if (e) return e;
    if (kp_is_one(&b->d)) return create_mul(kfrac_expr(b, &F->t1), ln);
    KFrac num = { kp_copy(&b->n), kp_one() };
    Expr* top = kp_is_one(&num.n) ? ln : create_mul(kfrac_expr(&num, &F->t1), ln);
    kfrac_free(&num);
    return create_div(top, lpoly_expr(&b->d, &F->t1));
}

// Σ_{k>=1} b_k·ln(u)^k (plus b_0 when with0), top power first
static Expr* poly2_expr(const K2Poly* B, const Form2* F, bool with0) {
    Expr* acc = create_const(0);
    for (int k = B->n - 1; k >= 1; k--)
        if (!kf_is_zero(&B->a[k])) acc = plus(acc, level_term(&B->a[k], F, k));
    if (with0 && B->n && !kf_is_zero(&B->a[0])) acc = plus(acc, kfrac_expr(&B->a[0], &F->t1));
    return acc;
}

// ∫ g0 over K1, the one-tower integrator
static RischResult integrate_k1(const KFrac* g0, const Tower* t) {
    ExpForm F0 = { *t, kp_copy(&g0->n), kp_copy(&g0->d), NULL };  // shares t
    Expr* e = kfrac_expr(g0, t);
    RischResult r = kp_is_one(&g0->d) ? integrate_poly(&F0, e) : integrate_fraction(&F0, e);
    free_expr(e);
    lpoly_free(&F0.f); lpoly_free(&F0.den);
    return r;
}

// ∫ P for a polynomial P in θ2: the levels, then ∫ g0 over K1 (f: the integrand, for a proof)
static RischResult integrate_poly2(const PolyIn* in, const Expr* f) {
    Levels L;
    LimStatus st = solve_levels(in, &L);
    RischResult r;
    if (st == LIM_NONE) r = (RischResult){ RISCH_NONELEMENTARY, expr_copy(f), LEVEL_FAILS };
    else if (st == LIM_UNSUPPORTED) r = failure(RISCH_UNSUPPORTED, "internal check failed: limited integration over the inner tower");
    else {
        K2Poly dB = k2_derivative(&L.B, in->E), G = k2_const(kf_copy(&L.g0), L.B.kind), s = k2_add(&dB, &G, 1), d = k2_add(&s, in->P, -1);
        bool same = d.n == 0;                              // the exact check: D(B) + g0 = P
        k2_free(&dB); k2_free(&G); k2_free(&s); k2_free(&d);
        r = same ? integrate_k1(&L.g0, &in->F->t1) : failure(RISCH_UNSUPPORTED, "internal check failed: the outer ln part does not differentiate back");
        if (r.status == RISCH_OK) r.result = plus(poly2_expr(&L.B, in->F, false), r.result);
    }
    k2_free(&L.B);
    kfrac_free(&L.g0);
    return r;
}

// --- denominators in θ2: Hermite reduction and the residue test over K1 ---
//
// Every factor of a denominator in θ2 = ln(u) is normal (θ2 is a primitive monomial, η ∉ D(K1)),
// so, as over Q(x) (§11.6): Bezout splits f = P + b/den, Hermite reduction with the derivation
// gives ∫ b/den = g + ∫ h/d (d square-free), the residues of h/d must be constants (their power
// sums Tr(s^k), s = h·D(d)⁻¹ mod d, all in Q), each rational residue c gives c·ln(gcd(d, h - c·D(d))),
// and f - D(g) - D(logs) is a polynomial in θ2 (checked exactly) for the levels above.

// D(n/d) = (D(n)·d - n·D(d))/d²
static K2Frac k2frac_derivative(const K2Frac* a, const Eta* E) {
    K2Poly dn = k2_derivative(&a->n, E), dd = k2_derivative(&a->d, E);
    K2Poly l = k2_mul(&dn, &a->d), r = k2_mul(&a->n, &dd), n = k2_add(&l, &r, -1);
    k2_free(&dn); k2_free(&dd); k2_free(&l); k2_free(&r);
    return k2frac_make(n, k2_mul(&a->d, &a->d));
}

static K2Poly k2_pow(const K2Poly* a, int k) {
    K2Poly r = k2_one(a->kind);
    for (int i = 0; i < k; i++) { K2Poly m = k2_mul(&r, a); k2_free(&r); r = m; }
    return r;
}

// c·a for a rational c
static K2Poly k2_scaled(const K2Poly* a, const Rational* c) {
    K2Poly r = k2_zeros(0, a->kind);
    for (int i = 0; i < a->n; i++) if (!kf_is_zero(&a->a[i])) k2_accumulate(&r, i, kf_scaled(&a->a[i], c, a->kind));
    return r;
}

// ∫ b/den = g + ∫ h/d
typedef struct {
    K2Frac g;
    K2Poly h, d;
} Reduced2;

// one Hermite step for V^(j+1) in U·V^(j+1): A/(U·V^(j+1)) = D(B/V^j) + (-j·C - U·D(B))/(U·V^j),
// with B·U·D(V) + C·V = -A/j
typedef struct {
    const K2Poly *U, *V, *UDV;
    int j;
    const Eta* E;
} Step2;

static void hermite2_step(Reduced2* R, K2Poly* A, const Step2* st) {
    const Eta* E = st->E;
    Rational mj = rat_from_si(-1, st->j), j = rat_from_int(-st->j);
    K2Poly rhs = k2_scaled(A, &mj);
    K2Bezout z = k2_bezout(st->UDV, st->V, &rhs);
    K2Frac term = k2frac_make(k2_copy(&z.s), k2_pow(st->V, st->j)), g = k2frac_add(&R->g, &term, 1);
    k2frac_free(&R->g);
    R->g = g;
    K2Poly DB = k2_derivative(&z.s, E), UDB = k2_mul(st->U, &DB), jC = k2_scaled(&z.t, &j);
    k2_free(A);
    *A = k2_add(&jC, &UDB, -1);
    k2_free(&rhs); k2bezout_free(&z); k2frac_free(&term);
    k2_free(&DB); k2_free(&UDB); k2_free(&jC);
    rat_clear(&mj); rat_clear(&j);
}

// Bronstein's quadratic Hermite reduction over K1[θ2], den monic
static Reduced2 hermite2(const K2Poly* b, const K2Poly* den, const Eta* E) {
    TowerKind kind = den->kind;
    Reduced2 R = { k2frac_zero(kind), k2_zeros(0, kind), k2_zeros(0, kind) };
    K2SquareFree sf = k2_square_free(den);
    K2Poly A = k2_copy(b), D = k2_copy(den);
    for (int i = 2; i <= sf.m; i++) {
        const K2Poly* V = &sf.f[i - 1];
        if (k2_deg(V) <= 0) continue;
        K2Poly Vi = k2_pow(V, i), U = k2_exact_div(&D, &Vi), DV = k2_derivative(V, E), UDV = k2_mul(&U, &DV);
        for (int j = i - 1; j >= 1; j--) hermite2_step(&R, &A, &(Step2){ &U, V, &UDV, j, E });
        k2_free(&D);
        D = k2_mul(&U, V);                                // V^i -> V
        k2_free(&Vi); k2_free(&U); k2_free(&DV); k2_free(&UDV);
    }
    K2Div q = k2_divmod(&A, &D);
    k2_free(&R.h); k2_free(&R.d);
    R.h = q.r;
    R.d = D;
    k2_free(&q.q); k2_free(&A);
    k2sf_free(&sf);
    return R;
}

static void reduced2_free(Reduced2* R) { k2frac_free(&R->g); k2_free(&R->h); k2_free(&R->d); }

// the integrand's form, η, and the Hermite reduction of its proper part
typedef struct {
    const Form2* F;
    const Eta* E;
    const Reduced2* R;
} Frac2;

// the coefficient of θ2^k (0 outside)
static KFrac k2_coeff(const K2Poly* p, int k) { return k >= 0 && k < p->n ? kf_copy(&p->a[k]) : kf_zero(); }

// P[j] = Σ θ_i^j over the roots of monic d, j < deg d (Newton's identities over K1)
static KFrac* root_power_sums2(const K2Poly* d) {
    TowerKind kind = d->kind;
    int n = k2_deg(d);
    KFrac* P = malloc((size_t)n * sizeof(KFrac));
    Rational nn = rat_from_int(n);
    QFn nq = qfn_rat(&nn);
    P[0] = kf_of_qfn(&nq, kind);
    qfn_free(&nq); rat_clear(&nn);
    for (int j = 1; j < n; j++) {                         // P_j = -(j·c_{n-j} + Σ_{i<j} c_{n-i}·P_{j-i})
        Rational jj = rat_from_int(j), m1 = rat_from_int(-1);
        KFrac cj = k2_coeff(d, n - j), acc = kf_scaled(&cj, &jj, kind);
        for (int i = 1; i < j; i++) {
            KFrac ci = k2_coeff(d, n - i), t = kfrac_mul(&ci, &P[j - i], kind), s2 = kfrac_add(&acc, &t, kind);
            kfrac_free(&acc); kfrac_free(&ci); kfrac_free(&t);
            acc = s2;
        }
        P[j] = kf_scaled(&acc, &m1, kind);
        kfrac_free(&cj); kfrac_free(&acc);
        rat_clear(&jj); rat_clear(&m1);
    }
    return P;
}

// Tr(g) = Σ_j g_j·P_j for deg g < deg d
static KFrac trace2(const K2Poly* g, const KFrac* P, int n) {
    KFrac acc = kf_zero();
    for (int j = 0; j < g->n && j < n; j++) {
        KFrac t = kfrac_mul(&g->a[j], &P[j], g->kind), s2 = kfrac_add(&acc, &t, g->kind);
        kfrac_free(&acc); kfrac_free(&t);
        acc = s2;
    }
    return acc;
}

typedef enum { RES_OK, RES_NONCONSTANT, RES_UNSUPPORTED } ResStatus;

// the power sums p_k = Tr(s^k), k = 1 .. n, of the residues; RES_NONCONSTANT when one is not in Q
static ResStatus residue_power_sums(const Reduced2* R, const Eta* E, Rational* p) {
    int n = k2_deg(&R->d);
    K2Poly Dd = k2_derivative(&R->d, E), one = k2_one(R->d.kind);
    K2Bezout z = k2_bezout(&Dd, &R->d, &one);             // z.s = D(d)⁻¹ mod d
    K2Poly hs = k2_mul(&R->h, &z.s);
    K2Div sd = k2_divmod(&hs, &R->d);                     // s = h·D(d)⁻¹ mod d
    KFrac* P = root_power_sums2(&R->d);
    K2Poly pw = k2_one(R->d.kind);
    ResStatus st = RES_OK;
    for (int k = 1; k <= n && st == RES_OK; k++) {
        K2Poly m = k2_mul(&pw, &sd.r);
        K2Div md = k2_divmod(&m, &R->d);
        k2_free(&pw); k2_free(&m); k2_free(&md.q);
        pw = md.r;
        KFrac tr = trace2(&pw, P, n);
        if (!kf_is_zero(&tr) && !kf_constant(&tr, &p[k - 1])) st = RES_NONCONSTANT;
        if (kf_is_zero(&tr)) rat_set_si(&p[k - 1], 0, 1);
        kfrac_free(&tr);
    }
    for (int j = 0; j < n; j++) kfrac_free(&P[j]);
    free(P);
    k2_free(&Dd); k2_free(&one); k2bezout_free(&z); k2_free(&hs); k2div_free(&sd); k2_free(&pw);
    return st;
}

// R(z) = Π (z - c_i) from the power sums p_1..p_n (Newton: k·e_k = Σ_{i=1}^{k} (-1)^(i-1)·e_{k-i}·p_i)
static Polynomial* residue_poly(const Rational* p, int n) {
    Rational* e = malloc((size_t)(n + 1) * sizeof(Rational));
    for (int k = 0; k <= n; k++) rat_init(&e[k]);
    rat_set_si(&e[0], 1, 1);
    Rational t, kk;
    rat_init(&t); rat_init(&kk);
    for (int k = 1; k <= n; k++) {
        for (int i = 1; i <= k; i++) {
            rat_mul(&t, &e[k - i], &p[i - 1]);
            if (i % 2) rat_add(&e[k], &e[k], &t); else rat_sub(&e[k], &e[k], &t);
        }
        rat_set_si(&kk, k, 1);
        rat_div(&e[k], &e[k], &kk);
    }
    Polynomial* R = create_polynomial();
    for (int k = 0; k <= n; k++) {
        if (k % 2) rat_neg(&t, &e[k]); else rat_set(&t, &e[k]);
        if (!rat_is_zero(&t)) add_term(R, &t, n - k);
    }
    for (int k = 0; k <= n; k++) rat_clear(&e[k]);
    free(e);
    rat_clear(&t); rat_clear(&kk);
    return R;
}

// the distinct roots of R, all rational, or false (a factor of degree >= 2)
typedef struct {
    Rational* c;
    int n;
} Roots;

static bool rational_roots(const Polynomial* R, Roots* out) {
    Factorization* fz = factorize(R);
    out->c = malloc((size_t)(fz->count + 1) * sizeof(Rational));
    out->n = 0;
    bool ok = true;
    for (int i = 0; i < fz->count && ok; i++) {
        Polynomial* q = fz->factors[i];
        ok = polynomial_degree(q) == 1;
        if (!ok) break;
        Rational a, b;
        rat_init(&a); rat_init(&b);
        poly_coeff(&a, q, 1);
        poly_coeff(&b, q, 0);
        rat_div(&b, &b, &a);
        rat_init(&out->c[out->n]);
        rat_neg(&out->c[out->n++], &b);
        rat_clear(&a); rat_clear(&b);
    }
    free_factorization(fz);
    return ok;
}

static void roots_free(Roots* r) { for (int i = 0; i < r->n; i++) rat_clear(&r->c[i]); free(r->c); }

// Σ c·ln(G_c), G_c = gcd(d, h - c·D(d)): the expression and its exact derivative
typedef struct {
    Expr* e;
    K2Frac deriv;
} LogPart2;

static void log_part2(const Frac2* X, const Roots* Z, LogPart2* out) {
    const Reduced2* R = X->R;
    const Form2* F = X->F;
    const Eta* E = X->E;
    K2Poly Dd = k2_derivative(&R->d, E);
    if (out->e) free_expr(out->e);
    out->e = create_const(0);
    k2frac_free(&out->deriv);
    out->deriv = k2frac_zero(R->d.kind);
    for (int i = 0; i < Z->n; i++) {
        K2Poly cD = k2_scaled(&Dd, &Z->c[i]), hc = k2_add(&R->h, &cD, -1), G = k2_gcd(&R->d, &hc);
        K2Poly DG = k2_derivative(&G, E), cDG = k2_scaled(&DG, &Z->c[i]);
        K2Frac t = k2frac_make(cDG, k2_copy(&G)), s2 = k2frac_add(&out->deriv, &t, 1);
        k2frac_free(&out->deriv); k2frac_free(&t);
        out->deriv = s2;
        Expr* term = create_mul(create_const_rat(&Z->c[i]), create_ln(poly2_expr(&G, F, true)));
        out->e = plus(out->e, term);
        k2_free(&cD); k2_free(&hc); k2_free(&G); k2_free(&DG);
    }
    k2_free(&Dd);
}

#define RESIDUES2_NOT_CONSTANT "has no elementary antiderivative (Risch: residues over the outer ln are not constants)"
#define FRAC2_LIMIT 6 // deg of a denominator in θ2 beyond which Euclid over K1 is not attempted

// g as an expression: num/den in ln(u)
static Expr* frac2_expr(const K2Frac* g, const Form2* F) {
    Expr* n = poly2_expr(&g->n, F, true);
    if (k2_deg(&g->d) <= 0) return n;
    return create_div(n, poly2_expr(&g->d, F, true));
}

// the log part of h/d with rational residues, or why not (out unchanged then)
static RischResult residues2(const Frac2* X, LogPart2* out) {
    const Reduced2* R = X->R;
    int n = k2_deg(&R->d);
    Rational* p = malloc((size_t)(n > 0 ? n : 1) * sizeof(Rational));
    for (int k = 0; k < n; k++) rat_init(&p[k]);
    ResStatus st = residue_power_sums(R, X->E, p);
    RischResult r = { RISCH_OK, NULL, NULL };
    if (st == RES_NONCONSTANT) r = (RischResult){ RISCH_NONELEMENTARY, NULL, RESIDUES2_NOT_CONSTANT };
    else {
        Polynomial* Rz = residue_poly(p, n);
        Roots Z;
        if (!rational_roots(Rz, &Z)) r = failure(RISCH_UNSUPPORTED, "residues over the outer ln that are irrational are not supported yet");
        else log_part2(X, &Z, out);
        roots_free(&Z);
        free_polynomial(Rz);
    }
    for (int k = 0; k < n; k++) rat_clear(&p[k]);
    free(p);
    return r;
}

// f - D(g) - D(logs) as a polynomial in θ2, or false
static bool remainder2(const Frac2* X, const K2Frac* logs, K2Poly* out) {
    const Form2* F = X->F;
    const Eta* E = X->E;
    const Reduced2* R = X->R;
    K2Frac Dg = k2frac_derivative(&R->g, E), a = k2frac_add(&F->f, &Dg, -1), b = k2frac_add(&a, logs, -1);
    bool ok = k2_deg(&b.d) == 0;
    if (ok) {
        KFrac ic = kfrac_inv(&b.d.a[0], b.d.kind);
        K2Poly icp = k2_const(ic, b.d.kind);
        *out = k2_mul(&b.n, &icp);
        k2_free(&icp);
    }
    k2frac_free(&Dg); k2frac_free(&a); k2frac_free(&b);
    return ok;
}

// ∫ f for a fraction f = N/D in θ2 (D nonconstant)
static RischResult integrate_frac2(const Form2* F, const Eta* E, const Expr* f) {
    if (k2_deg(&F->f.d) > FRAC2_LIMIT) return failure(RISCH_UNSUPPORTED, "a denominator in the outer ln of high degree is not supported yet");
    K2Div qd = k2_divmod(&F->f.n, &F->f.d);              // f = Q + R/D
    Reduced2 R = hermite2(&qd.r, &F->f.d, E);
    LogPart2 L = { NULL, k2frac_zero(R.d.kind) };
    Frac2 X = { F, E, &R };
    RischResult r = k2_deg(&R.d) > 0 && R.h.n ? residues2(&X, &L) : (RischResult){ RISCH_OK, NULL, NULL };
    if (r.status == RISCH_NONELEMENTARY) r.result = frac2_expr(&(K2Frac){ R.h, R.d }, F);
    K2Poly P;
    if (r.status == RISCH_OK && !remainder2(&X, &L.deriv, &P)) r = failure(RISCH_UNSUPPORTED, "internal check failed: the remainder over the outer ln is not a polynomial");
    else if (r.status == RISCH_OK) {
        r = integrate_poly2(&(PolyIn){ F, E, &P }, f);
        if (r.status == RISCH_OK) {
            if (R.g.n.n) r.result = plus(r.result, frac2_expr(&R.g, F));
            if (L.e) { r.result = plus(r.result, L.e); L.e = NULL; }
        }
        k2_free(&P);
    }
    if (L.e) free_expr(L.e);
    k2frac_free(&L.deriv);
    reduced2_free(&R);
    k2div_free(&qd);
    return r;
}

RischResult integrate_tower2(const Expr* f) {
    Form2 F;
    if (!form2_of(f, &F)) return failure(RISCH_UNSUPPORTED, F.why);
    Eta E;
    RischResult r;
    if (!eta_init(&F.u, &F.t1, &E)) r = failure(RISCH_UNSUPPORTED, "internal check failed: the outer ln's derivative over the inner tower");
    else if (eta_in_dk1(&E)) r = failure(RISCH_UNSUPPORTED, "the outer ln lies in the inner tower (such as ln(exp(x)) or ln(2x) beside ln(x)); this is not supported yet");
    else if (k2_deg(&F.f.d) == 0) {                       // a polynomial in θ2 (the denominator is monic: 1)
        r = integrate_poly2(&(PolyIn){ &F, &E, &F.f.n }, f);
    } else r = integrate_frac2(&F, &E, f);
    eta_free(&E);
    form2_free(&F);
    return r;
}
