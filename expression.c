#include "expression.h"
#include "polynomial.h"
#include "factor.h"
#include "strbuf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "roots.h"
#include "qsurd.h"

// --- constructors ---

Expr* create_const_rat(const Rational* value) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_CONST;
    rat_init_set(&e->constant, value);
    return e;
}

Expr* create_const(long value) {
    Rational r = rat_from_int(value);
    Expr* e = create_const_rat(&r);
    rat_clear(&r);
    return e;
}

Expr* create_const_str(const char* literal) {
    Rational r;
    rat_init(&r);
    if (!rat_parse_literal(&r, literal)) mpq_set_si(r.val, 0, 1);
    Expr* e = create_const_rat(&r);
    rat_clear(&r);
    return e;
}

Expr* create_var(const char* name) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_VAR;
    e->var_name = strdup(name);
    return e;
}

static Expr* create_binop(BinOp op, Expr* left, Expr* right) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_BINOP;
    e->binop.op = op;
    e->binop.left = left;
    e->binop.right = right;
    return e;
}

Expr* create_add(Expr* left, Expr* right) { return create_binop(OP_ADD, left, right); }
Expr* create_sub(Expr* left, Expr* right) { return create_binop(OP_SUB, left, right); }
Expr* create_mul(Expr* left, Expr* right) { return create_binop(OP_MUL, left, right); }
Expr* create_div(Expr* numerator, Expr* denominator) { return create_binop(OP_DIV, numerator, denominator); }
Expr* create_pow(Expr* base, Expr* exponent) { return create_binop(OP_POW, base, exponent); }

Expr* create_neg(Expr* operand) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_NEG;
    e->operand = operand;
    return e;
}

Expr* create_poly_expr(Polynomial* poly) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_POLY;
    e->poly = poly;
    return e;
}

Expr* create_rational_fn(Polynomial* num, Polynomial* den) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_RATIONAL_FN;
    e->ratfn.num = num;
    e->ratfn.den = den;
    return e;
}

Expr* create_ln(Expr* arg) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_LN;
    e->ln_arg = arg;
    return e;
}

static Expr* create_unary_fn(ExprType type, Expr* arg) {
    Expr* e = malloc(sizeof(Expr));
    e->type = type;
    e->fn_arg = arg;
    return e;
}

Expr* create_atan(Expr* arg) { return create_unary_fn(EXPR_ATAN, arg); }
Expr* create_exp(Expr* arg) { return create_unary_fn(EXPR_EXP, arg); }
Expr* create_sin(Expr* arg) { return create_unary_fn(EXPR_SIN, arg); }
Expr* create_cos(Expr* arg) { return create_unary_fn(EXPR_COS, arg); }

Expr* create_rootsum(Polynomial* q, Polynomial** s, int s_deg) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_ROOTSUM;
    e->rootsum.q = q;
    e->rootsum.s = s;
    e->rootsum.s_deg = s_deg;
    return e;
}
Expr* create_sqrt(Expr* arg) { return create_unary_fn(EXPR_SQRT, arg); }

Expr* create_trig(TrigKind kind, const Rational* angle) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_TRIG;
    e->trig.kind = kind;
    rat_init_set(&e->trig.angle, angle);
    return e;
}

Expr* create_pi(void) {
    Rational one = rat_one_val();
    Expr* e = create_trig(TRIG_NONE, &one);
    rat_clear(&one);
    return e;
}

Expr* expr_copy(const Expr* e) {
    switch (e->type) {
        case EXPR_CONST: return create_const_rat(&e->constant);
        case EXPR_VAR: return create_var(e->var_name);
        case EXPR_BINOP: return create_binop(e->binop.op, expr_copy(e->binop.left), expr_copy(e->binop.right));
        case EXPR_NEG: return create_neg(expr_copy(e->operand));
        case EXPR_POLY: return create_poly_expr(poly_copy(e->poly));
        case EXPR_RATIONAL_FN: return create_rational_fn(poly_copy(e->ratfn.num), poly_copy(e->ratfn.den));
        case EXPR_LN: return create_ln(expr_copy(e->ln_arg));
        case EXPR_ATAN: return create_atan(expr_copy(e->fn_arg));
        case EXPR_EXP: return create_exp(expr_copy(e->fn_arg));
        case EXPR_SIN: return create_sin(expr_copy(e->fn_arg));
        case EXPR_COS: return create_cos(expr_copy(e->fn_arg));
        case EXPR_SQRT: return create_sqrt(expr_copy(e->fn_arg));
        case EXPR_TRIG: return create_trig(e->trig.kind, &e->trig.angle);
        case EXPR_ROOTSUM: {
            Polynomial** s = malloc((size_t)(e->rootsum.s_deg + 1) * sizeof(Polynomial*));
            for (int k = 0; k <= e->rootsum.s_deg; k++) s[k] = poly_copy(e->rootsum.s[k]);
            return create_rootsum(poly_copy(e->rootsum.q), s, e->rootsum.s_deg);
        }
    }
    return NULL;
}

// a node of this type anywhere in e
static bool has_type(const Expr* e, ExprType t) {
    if (e->type == t) return true;
    switch (e->type) {
        case EXPR_BINOP: return has_type(e->binop.left, t) || has_type(e->binop.right, t);
        case EXPR_NEG: return has_type(e->operand, t);
        case EXPR_LN: return has_type(e->ln_arg, t);
        case EXPR_ATAN: case EXPR_SQRT: case EXPR_EXP: case EXPR_SIN: case EXPR_COS: return has_type(e->fn_arg, t);
        default: return false;
    }
}

bool expr_has_exp(const Expr* e) { return has_type(e, EXPR_EXP); }
bool expr_has_ln(const Expr* e) { return has_type(e, EXPR_LN); }
bool expr_has_trig(const Expr* e) { return has_type(e, EXPR_SIN) || has_type(e, EXPR_COS); }

bool expr_has_rootsum(const Expr* e) {
    if (e->type == EXPR_ROOTSUM) return true;
    if (e->type == EXPR_BINOP) return expr_has_rootsum(e->binop.left) || expr_has_rootsum(e->binop.right);
    return e->type == EXPR_NEG && expr_has_rootsum(e->operand);
}

// --- memory ---

void free_expr(Expr* expr) {
    if (!expr) return;
    switch (expr->type) {
        case EXPR_VAR:
            free(expr->var_name);
            break;
        case EXPR_BINOP:
            free_expr(expr->binop.left);
            free_expr(expr->binop.right);
            break;
        case EXPR_NEG:
            free_expr(expr->operand);
            break;
        case EXPR_POLY:
            free_polynomial(expr->poly);
            break;
        case EXPR_RATIONAL_FN:
            free_polynomial(expr->ratfn.num);
            free_polynomial(expr->ratfn.den);
            break;
        case EXPR_LN:
            free_expr(expr->ln_arg);
            break;
        case EXPR_ATAN:
        case EXPR_SQRT:
        case EXPR_EXP:
        case EXPR_SIN:
        case EXPR_COS:
            free_expr(expr->fn_arg);
            break;
        case EXPR_ROOTSUM:
            free_polynomial(expr->rootsum.q);
            for (int k = 0; k <= expr->rootsum.s_deg; k++) free_polynomial(expr->rootsum.s[k]);
            free(expr->rootsum.s);
            break;
        case EXPR_TRIG:
            rat_clear(&expr->trig.angle);
            break;
        case EXPR_CONST:
            rat_clear(&expr->constant);
            break;
    }
    free(expr);
}

// --- output ---

// check if an expression is a negative constant
static bool is_neg_const(const Expr* e) {
    return e->type == EXPR_CONST && rat_is_negative(&e->constant);
}

// check if expr is the constant 1
static bool is_one_const(const Expr* e) {
    return e->type == EXPR_CONST && rat_is_one(&e->constant);
}

// check if expr is the constant -1
static bool is_neg_one_const(const Expr* e) {
    if (e->type != EXPR_CONST) return false;
    Rational neg1 = rat_from_int(-1);
    bool result = rat_eq(&e->constant, &neg1);
    rat_clear(&neg1);
    return result;
}

// check if an expression is a negative-leading term (neg const, or neg const * something)
static bool is_negative_term(const Expr* e) {
    if (is_neg_const(e)) return true;
    if (e->type == EXPR_NEG) return true;
    if (e->type == EXPR_BINOP && e->binop.op == OP_MUL && is_neg_const(e->binop.left)) return true;
    return false;
}

// a positive rational, or a positive element p + q·sqrt(m) (decided exactly)
static bool positive_constant(const Expr* e) {
    if (e->type == EXPR_CONST) return mpq_sgn(e->constant.val) > 0;
    QSurd a;
    qs_init(&a);
    bool pos = qs_parse(e, &a) && qs_sign(&a) > 0;
    qs_clear(&a);
    return pos;
}

// positive for every real x by its shape: positive constants and exp(...) combined by + · /
// (exp(x) + 1, 2·exp(x)), so ln of it needs no abs
static bool evidently_positive(const Expr* e) {
    switch (e->type) {
        case EXPR_EXP: return true;
        case EXPR_POLY: return polynomial_degree(e->poly) == 0 && !rat_is_negative(&e->poly->head->coeff);
        case EXPR_BINOP:
            return (e->binop.op == OP_ADD || e->binop.op == OP_MUL || e->binop.op == OP_DIV)
                && evidently_positive(e->binop.left) && evidently_positive(e->binop.right);
        default: return positive_constant(e);
    }
}

static void expr_to_string_impl(const Expr* expr, StrBuf* sb);

// emit a term with its sign flipped: -1/2 * ln(abs(x)) → 1/2 * ln(abs(x))
static void emit_negated_term(const Expr* expr, StrBuf* sb) {
    if (is_neg_const(expr)) {
        Rational pos; rat_init(&pos);
        rat_neg(&pos, &expr->constant);
        char* s = rat_to_string(&pos);
        sb_append(sb, s);
        free(s);
        rat_clear(&pos);
        return;
    }
    if (expr->type == EXPR_NEG) {
        expr_to_string_impl(expr->operand, sb);
        return;
    }
    if (expr->type == EXPR_BINOP && expr->binop.op == OP_MUL
        && is_neg_const(expr->binop.left)) {
        Rational pos; rat_init(&pos);
        rat_neg(&pos, &expr->binop.left->constant);
        if (rat_is_one(&pos)) {
            expr_to_string_impl(expr->binop.right, sb);
        } else {
            char* s = rat_to_string(&pos);
            sb_append(sb, s);
            free(s);
            sb_append(sb, " * ");
            expr_to_string_impl(expr->binop.right, sb);
        }
        rat_clear(&pos);
        return;
    }
    expr_to_string_impl(expr, sb);
}

// a polynomial whose leading coefficient is negative (it prints with a leading "-")
static bool leading_minus_poly(const Expr* e) {
    return e->type == EXPR_POLY && e->poly->head && rat_is_negative(&e->poly->head->coeff);
}

// a product led by a one-term polynomial with a negative coefficient: -1/3x^3 * ln(abs(x))
static bool negative_monomial_product(const Expr* e) {
    return e->type == EXPR_BINOP && (e->binop.op == OP_MUL || e->binop.op == OP_DIV)
        && leading_minus_poly(e->binop.left) && !e->binop.left->poly->head->next;
}

// emit a + b, absorbing negative right-hand sides as subtraction
static void emit_add(const Expr* left, const Expr* right, StrBuf* sb) {
    expr_to_string_impl(left, sb);
    if (leading_minus_poly(right) || negative_monomial_product(right)) {   // a + -x^2 + 3 -> a - x^2 + 3
        StrBuf rs = sb_new();
        expr_to_string_impl(right, &rs);
        char* str = sb_finish(&rs);
        sb_append(sb, " - ");
        sb_append(sb, str + 1);                           // drop the leading '-'
        free(str);
    } else if (is_negative_term(right)) {
        sb_append(sb, " - ");
        emit_negated_term(right, sb);
    } else {
        sb_append(sb, " + ");
        expr_to_string_impl(right, sb);
    }
}

// true for a sum that needs parentheses as a factor: a+b, a-b, or a multi-term polynomial
static bool is_sum(const Expr* e) {
    if (e->type == EXPR_BINOP) return e->binop.op == OP_ADD || e->binop.op == OP_SUB;
    return e->type == EXPR_POLY && e->poly->head && e->poly->head->next;
}

// emit an operand, parenthesized when it is a sum
static void emit_factor(const Expr* e, StrBuf* sb) {
    bool p = is_sum(e);
    if (p) sb_append(sb, "(");
    expr_to_string_impl(e, sb);
    if (p) sb_append(sb, ")");
}

// emit a * b, eliding coefficients of 1 and -1
static void emit_mul(const Expr* left, const Expr* right, StrBuf* sb) {
    if (is_one_const(left)) {
        emit_factor(right, sb);
    } else if (is_neg_one_const(left)) {
        sb_append(sb, "-");
        emit_factor(right, sb);
    } else {
        emit_factor(left, sb);
        sb_append(sb, " * ");
        emit_factor(right, sb);
    }
}

// sin(u)/cos(u), printed as tan(u)
static bool is_tan(const Expr* e) {
    if (e->type != EXPR_BINOP || e->binop.op != OP_DIV) return false;
    const Expr *l = e->binop.left, *r = e->binop.right;
    if (l->type != EXPR_SIN || r->type != EXPR_COS) return false;
    char *a = expr_to_string(l->fn_arg), *b = expr_to_string(r->fn_arg);
    bool same = strcmp(a, b) == 0;
    free(a); free(b);
    return same;
}

// emit a / b with minimal parenthesization
static void emit_div(const Expr* left, const Expr* right, StrBuf* sb) {
    bool lp = is_sum(left);
    bool rp = (right->type == EXPR_BINOP && !is_tan(right)) || is_sum(right);
    if (lp) sb_append(sb, "(");
    expr_to_string_impl(left, sb);
    if (lp) sb_append(sb, ")");
    sb_append(sb, " / ");
    if (rp) sb_append(sb, "(");
    expr_to_string_impl(right, sb);
    if (rp) sb_append(sb, ")");
}

// emit base^exp, rendering negative exponents as fractions
static void emit_pow(const Expr* base, const Expr* exp, StrBuf* sb) {
    if (is_neg_const(exp)) {
        Rational pos_exp; rat_init(&pos_exp);
        rat_neg(&pos_exp, &exp->constant);
        sb_append(sb, "1/(");
        expr_to_string_impl(base, sb);
        sb_append(sb, ")");
        if (!rat_is_one(&pos_exp)) {
            sb_append(sb, "^");
            char* s = rat_to_string(&pos_exp);
            sb_append(sb, s);
            free(s);
        }
        rat_clear(&pos_exp);
        return;
    }
    bool needs_parens = (base->type == EXPR_BINOP && !is_tan(base)) || base->type == EXPR_POLY;
    if (needs_parens) sb_append(sb, "(");
    expr_to_string_impl(base, sb);
    if (needs_parens) sb_append(sb, ")");
    sb_append(sb, "^");
    bool frac = exp->type == EXPR_CONST && mpz_cmp_ui(mpq_denref(exp->constant.val), 1) != 0; // 2^(1/6), not 2^1/6
    if (frac) sb_append(sb, "(");
    expr_to_string_impl(exp, sb);
    if (frac) sb_append(sb, ")");
}

// cos(pi/7), sin(3*pi/10), ...
static void emit_trig(const Expr* e, StrBuf* sb) {
    static const char* names[] = { "", "cos", "sin", "cot", "csc" };
    if (e->trig.kind == TRIG_NONE) { sb_append(sb, "pi"); return; } // create_pi: angle 1
    const mpq_t* a = &e->trig.angle.val;
    sb_appendf(sb, "%s(", names[e->trig.kind]);
    if (mpz_cmp_ui(mpq_numref(*a), 1) != 0) { sb_append_owned(sb, mpz_get_str(NULL, 10, mpq_numref(*a))); sb_append(sb, "*"); }
    sb_append(sb, "pi");
    if (mpz_cmp_ui(mpq_denref(*a), 1) != 0) { sb_append(sb, "/"); sb_append_owned(sb, mpz_get_str(NULL, 10, mpq_denref(*a))); }
    sb_append(sb, ")");
}

// one x^k term of a bivariate S(t, x) whose coefficient c is a polynomial in t
static void emit_tx_term(StrBuf* sb, const Polynomial* c, int k) {
    bool first = sb->len == 0 || sb->buf[sb->len - 1] == '(';
    bool mono = c->head && !c->head->next;          // single term: fold its sign
    bool neg = mono && rat_is_negative(&c->head->coeff);
    if (!first) sb_append(sb, neg ? " - " : " + ");
    else if (neg) sb_append(sb, "-");
    Polynomial* shown = mono ? create_polynomial() : poly_copy(c);
    if (mono) { Rational a; rat_init(&a); rat_abs(&a, &c->head->coeff); add_term(shown, &a, c->head->exponent); rat_clear(&a); }
    bool unit = mono && c->head->exponent == 0 && rat_is_one(&shown->head->coeff);
    if (k == 0 || !unit) {
        if (!mono) sb_append(sb, "(");
        sb_append_owned(sb, poly_to_string_var(shown, "t"));
        if (!mono) sb_append(sb, ")");
        if (k > 0) sb_append(sb, "*");
    }
    if (k > 0) sb_append(sb, "x");
    if (k > 1) sb_appendf(sb, "^%d", k);
    free_polynomial(shown);
}

// RootSum(q(t), t -> t * ln(S(t, x)))
static void emit_rootsum(const Expr* e, StrBuf* sb) {
    sb_append(sb, "RootSum(");
    sb_append_owned(sb, poly_to_string_var(e->rootsum.q, "t"));
    sb_append(sb, ", t -> t * ln(");
    for (int k = e->rootsum.s_deg; k >= 0; k--)
        if (!poly_is_zero(e->rootsum.s[k])) emit_tx_term(sb, e->rootsum.s[k], k);
    sb_append(sb, "))");
}

static void expr_to_string_impl(const Expr* expr, StrBuf* sb) {
    if (!expr) { sb_append(sb, "NULL"); return; }

    switch (expr->type) {
        case EXPR_CONST: {
            char* s = rat_to_string(&expr->constant);
            sb_append(sb, s);
            free(s);
            break;
        }
        case EXPR_VAR:
            sb_append(sb, expr->var_name);
            break;
        case EXPR_NEG: {
            ExprType t = expr->operand->type;              // a function call needs no parentheses
            bool call = t == EXPR_EXP || t == EXPR_SIN || t == EXPR_COS || t == EXPR_LN || t == EXPR_ATAN || t == EXPR_SQRT || t == EXPR_TRIG
                     || (t == EXPR_BINOP && expr->operand->binop.op != OP_ADD && expr->operand->binop.op != OP_SUB);
            sb_append(sb, call ? "-" : "-(");
            expr_to_string_impl(expr->operand, sb);
            if (!call) sb_append(sb, ")");
            break;
        }
        case EXPR_BINOP:
            switch (expr->binop.op) {
                case OP_ADD: emit_add(expr->binop.left, expr->binop.right, sb); break;
                case OP_SUB:
                    expr_to_string_impl(expr->binop.left, sb);
                    sb_append(sb, " - ");
                    expr_to_string_impl(expr->binop.right, sb);
                    break;
                case OP_MUL: emit_mul(expr->binop.left, expr->binop.right, sb); break;
                case OP_DIV:
                    if (is_tan(expr)) {
                        sb_append(sb, "tan(");
                        expr_to_string_impl(expr->binop.left->fn_arg, sb);
                        sb_append(sb, ")");
                    } else emit_div(expr->binop.left, expr->binop.right, sb);
                    break;
                case OP_POW: emit_pow(expr->binop.left, expr->binop.right, sb); break;
            }
            break;
        case EXPR_POLY: {
            char* s = poly_to_string(expr->poly);
            sb_append(sb, s);
            free(s);
            break;
        }
        case EXPR_RATIONAL_FN: {
            sb_append(sb, "(");
            char* ns = poly_to_string(expr->ratfn.num);
            sb_append(sb, ns);
            free(ns);
            sb_append(sb, ") / (");
            char* ds = poly_to_string(expr->ratfn.den);
            sb_append(sb, ds);
            free(ds);
            sb_append(sb, ")");
            break;
        }
        case EXPR_LN:
            if (evidently_positive(expr->ln_arg)) {
                sb_append(sb, "ln(");                    // ln(3/2), ln(3 + 2 * sqrt(2)), ln(exp(x) + 1): no abs
                expr_to_string_impl(expr->ln_arg, sb);
                sb_append(sb, ")");
                break;
            }
            sb_append(sb, "ln(abs(");
            expr_to_string_impl(expr->ln_arg, sb);
            sb_append(sb, "))");
            break;
        case EXPR_ATAN:
        case EXPR_SQRT:
        case EXPR_EXP:
        case EXPR_SIN:
        case EXPR_COS: {
            static const char* names[] = { [EXPR_ATAN] = "atan(", [EXPR_SQRT] = "sqrt(", [EXPR_EXP] = "exp(", [EXPR_SIN] = "sin(", [EXPR_COS] = "cos(" };
            sb_append(sb, names[expr->type]);
            expr_to_string_impl(expr->fn_arg, sb);
            sb_append(sb, ")");
            break;
        }
        case EXPR_ROOTSUM:
            emit_rootsum(expr, sb);
            break;
        case EXPR_TRIG:
            emit_trig(expr, sb);
            break;
    }
}

char* expr_to_string(const Expr* expr) {
    StrBuf sb = sb_new();
    expr_to_string_impl(expr, &sb);
    return sb_finish(&sb);
}

void print_expr(const Expr* expr) {
    char* s = expr_to_string(expr);
    printf("%s", s);
    free(s);
}

// --- canonicalization ---

static Expr* wrap_poly(Polynomial* poly) {
    return create_poly_expr(poly);
}

static Expr* wrap_ratfn(Polynomial* num, Polynomial* den) {
    if (poly_is_zero(den)) { // division by zero is an error, not a value
        free_polynomial(num);
        free_polynomial(den);
        return NULL;
    }
    if (polynomial_degree(den) == 0) { // constant denominator: scale the numerator
        Rational inv;
        rat_init(&inv);
        mpq_inv(inv.val, den->head->coeff.val);
        Polynomial* c = poly_constant(&inv);
        Polynomial* scaled = poly_multiply(num, c);
        free_polynomial(c); free_polynomial(num); free_polynomial(den);
        rat_clear(&inv);
        return wrap_poly(scaled);
    }
    Polynomial* g = poly_gcd(num, den);
    if (!poly_is_zero(g) && polynomial_degree(g) > 0) {
        Polynomial* rn = poly_divide(num, g);
        Polynomial* rd = poly_divide(den, g);
        free_polynomial(num);
        free_polynomial(den);
        free_polynomial(g);
        num = rn;
        den = rd;
    } else {
        free_polynomial(g);
    }
    if (den->head != NULL && den->head->next == NULL &&
        den->head->exponent == 0 && rat_is_one(&den->head->coeff)) {
        free_polynomial(den);
        return wrap_poly(num);
    }
    return create_rational_fn(num, den);
}

// largest polynomial degree a power may produce (guards against x^999999999)
#define MAX_POW_DEGREE 100000

// exact integer exponent from a constant polynomial; false if not an integer or too large
static bool exact_int_exponent(const Polynomial* rp, long* out) {
    if (rp->head == NULL) { *out = 0; return true; }
    if (rp->head->exponent != 0 || rp->head->next != NULL) return false;
    const mpq_t* q = &rp->head->coeff.val;
    if (mpz_cmp_ui(mpq_denref(*q), 1) != 0) return false; // x^(1/2) is not a polynomial
    if (!mpz_fits_slong_p(mpq_numref(*q))) return false;
    *out = mpz_get_si(mpq_numref(*q));
    return true;
}

// lp^rp for a constant integer rp; consumes both
static Expr* canon_int_pow(Polynomial* lp, Polynomial* rp) {
    long exp;
    bool ok = exact_int_exponent(rp, &exp);
    free_polynomial(rp);
    long deg = polynomial_degree(lp) > 0 ? polynomial_degree(lp) : 1;
    long mag = exp < 0 ? -exp : exp;
    if (!ok || mag > MAX_POW_DEGREE / deg) { free_polynomial(lp); return NULL; }
    Polynomial* powered = poly_pow(lp, (unsigned)mag);
    free_polynomial(lp);
    if (exp >= 0) return wrap_poly(powered);
    Rational one = rat_one_val();
    Polynomial* num = poly_constant(&one);
    rat_clear(&one);
    return wrap_ratfn(num, powered);
}

static Expr* canon_poly_op(BinOp op, Polynomial* lp, Polynomial* rp) {
    switch (op) {
        case OP_ADD: {
            Polynomial* result = poly_add(lp, rp);
            free_polynomial(lp); free_polynomial(rp);
            return wrap_poly(result);
        }
        case OP_SUB: {
            Polynomial* result = poly_subtract(lp, rp);
            free_polynomial(lp); free_polynomial(rp);
            return wrap_poly(result);
        }
        case OP_MUL: {
            Polynomial* result = poly_multiply(lp, rp);
            free_polynomial(lp); free_polynomial(rp);
            return wrap_poly(result);
        }
        case OP_DIV:
            return wrap_ratfn(lp, rp);
        case OP_POW:
            return canon_int_pow(lp, rp);
    }
    free_polynomial(lp); free_polynomial(rp);
    return NULL;
}

static Polynomial* extract_poly(Expr* e) {
    if (e->type == EXPR_POLY) {
        Polynomial* p = e->poly;
        e->poly = NULL;
        return p;
    }
    return NULL;
}

static bool extract_ratfn(Expr* e, Polynomial** num, Polynomial** den) {
    if (e->type == EXPR_POLY) {
        *num = e->poly;
        e->poly = NULL;
        *den = create_polynomial();
        Rational one_r = rat_one_val();
        add_term(*den, &one_r, 0);
        rat_clear(&one_r);
        return true;
    }
    if (e->type == EXPR_RATIONAL_FN) {
        *num = e->ratfn.num;
        *den = e->ratfn.den;
        e->ratfn.num = NULL;
        e->ratfn.den = NULL;
        return true;
    }
    return false;
}

Expr* canonicalize(Expr* expr) {
    if (!expr) return NULL;

    switch (expr->type) {
        case EXPR_CONST: {
            Polynomial* p = create_polynomial();
            add_term(p, &expr->constant, 0);
            return wrap_poly(p);
        }
        case EXPR_VAR: {
            Polynomial* p = create_polynomial();
            Rational one_r = rat_one_val();
            add_term(p, &one_r, 1);
            rat_clear(&one_r);
            return wrap_poly(p);
        }
        case EXPR_NEG: {
            Expr* inner = canonicalize(expr->operand);
            if (!inner) return NULL;
            if (inner->type == EXPR_POLY) {
                Polynomial* neg = create_polynomial();
                Term* t = inner->poly->head;
                while (t) {
                    Rational neg_c;
                    rat_init(&neg_c);
                    rat_neg(&neg_c, &t->coeff);
                    add_term(neg, &neg_c, t->exponent);
                    rat_clear(&neg_c);
                    t = t->next;
                }
                free_expr(inner);
                return wrap_poly(neg);
            }
            if (inner->type == EXPR_RATIONAL_FN) {
                Polynomial* neg_num = create_polynomial();
                Term* t = inner->ratfn.num->head;
                while (t) {
                    Rational neg_c;
                    rat_init(&neg_c);
                    rat_neg(&neg_c, &t->coeff);
                    add_term(neg_num, &neg_c, t->exponent);
                    rat_clear(&neg_c);
                    t = t->next;
                }
                Polynomial* den = poly_copy(inner->ratfn.den);
                free_expr(inner);
                return wrap_ratfn(neg_num, den);
            }
            free_expr(inner);
            return NULL;
        }
        case EXPR_BINOP: {
            Expr* left = canonicalize(expr->binop.left);
            Expr* right = canonicalize(expr->binop.right);
            if (!left || !right) {
                free_expr(left); free_expr(right);
                return NULL;
            }
            if (left->type == EXPR_POLY && right->type == EXPR_POLY) {
                Polynomial* lp = extract_poly(left);
                Polynomial* rp = extract_poly(right);
                free(left); free(right);
                return canon_poly_op(expr->binop.op, lp, rp);
            }
            Polynomial *ln, *ld, *rn, *rd;
            if (!extract_ratfn(left, &ln, &ld) || !extract_ratfn(right, &rn, &rd)) {
                free_expr(left); free_expr(right);
                return NULL;
            }
            free(left); free(right);

            switch (expr->binop.op) {
                case OP_ADD: case OP_SUB: {
                    Polynomial* ad = poly_multiply(ln, rd);
                    Polynomial* cb = poly_multiply(rn, ld);
                    Polynomial* num = (expr->binop.op == OP_ADD) ?
                        poly_add(ad, cb) : poly_subtract(ad, cb);
                    Polynomial* den = poly_multiply(ld, rd);
                    free_polynomial(ad); free_polynomial(cb);
                    free_polynomial(ln); free_polynomial(ld);
                    free_polynomial(rn); free_polynomial(rd);
                    return wrap_ratfn(num, den);
                }
                case OP_MUL: {
                    Polynomial* num = poly_multiply(ln, rn);
                    Polynomial* den = poly_multiply(ld, rd);
                    free_polynomial(ln); free_polynomial(ld);
                    free_polynomial(rn); free_polynomial(rd);
                    return wrap_ratfn(num, den);
                }
                case OP_DIV: {
                    Polynomial* num = poly_multiply(ln, rd);
                    Polynomial* den = poly_multiply(ld, rn);
                    free_polynomial(ln); free_polynomial(ld);
                    free_polynomial(rn); free_polynomial(rd);
                    return wrap_ratfn(num, den);
                }
                case OP_POW: {
                    free_polynomial(rd);
                    Polynomial* result_num = canon_poly_op(OP_POW, ln, rn) ? NULL : NULL;
                    (void)result_num;
                    free_polynomial(ld);
                    return NULL;
                }
            }
            return NULL;
        }
        case EXPR_POLY:
            return wrap_poly(poly_copy(expr->poly));
        case EXPR_RATIONAL_FN:
            return wrap_ratfn(poly_copy(expr->ratfn.num), poly_copy(expr->ratfn.den));
        case EXPR_LN:
        case EXPR_ATAN:
        case EXPR_SQRT:
        case EXPR_ROOTSUM:
        case EXPR_TRIG:
        case EXPR_EXP:                                    // not a rational function (risch.c)
        case EXPR_SIN:
        case EXPR_COS:
            return NULL;
    }
    return NULL;
}

Expr* simplify(Expr* expr) {
    return canonicalize(expr);
}

Expr* factor(Expr* expr) {
    if (!expr) return NULL;
    Expr* canonical = canonicalize(expr);
    if (!canonical || canonical->type != EXPR_POLY) {
        free_expr(canonical);
        return NULL;
    }

    Factorization* f = factorize(canonical->poly);
    free_expr(canonical);
    if (!f) return NULL;

    Expr* result = NULL;

    if (f->count == 0) {
        result = create_const_rat(&f->content);
        free_factorization(f);
        return result;
    }

    for (int i = 0; i < f->count; i++) {
        Expr* factor_expr = create_poly_expr(poly_copy(f->factors[i]));
        if (f->multiplicities[i] > 1) {
            Rational mi = rat_from_int(f->multiplicities[i]);
            factor_expr = create_pow(factor_expr, create_const_rat(&mi));
            rat_clear(&mi);
        }
        if (result == NULL)
            result = factor_expr;
        else
            result = create_mul(result, factor_expr);
    }

    if (!rat_is_one(&f->content)) {
        result = create_mul(create_const_rat(&f->content), result);
    }

    free_factorization(f);
    return result;
}

Expr* differentiate(Expr* expr, const char* variable) {
    if (!expr) return NULL;
    (void)variable;

    Expr* canonical = canonicalize(expr);
    if (!canonical) return NULL;

    if (canonical->type == EXPR_POLY) {
        Polynomial* deriv = poly_derivative(canonical->poly);
        free_expr(canonical);
        return wrap_poly(deriv);
    }

    if (canonical->type == EXPR_RATIONAL_FN) {
        Polynomial* f = canonical->ratfn.num;
        Polynomial* g = canonical->ratfn.den;
        Polynomial* fp = poly_derivative(f);
        Polynomial* gp = poly_derivative(g);
        Polynomial* fp_g = poly_multiply(fp, g);
        Polynomial* f_gp = poly_multiply(f, gp);
        Polynomial* num = poly_subtract(fp_g, f_gp);
        Polynomial* den = poly_multiply(g, g);
        free_polynomial(fp); free_polynomial(gp);
        free_polynomial(fp_g); free_polynomial(f_gp);
        free_expr(canonical);
        return wrap_ratfn(num, den);
    }

    free_expr(canonical);
    return NULL;
}

// --- numeric evaluation ---

static double eval_binop(const Expr* e, double x, bool* ok) {
    double l = expr_eval(e->binop.left, x, ok), r = expr_eval(e->binop.right, x, ok);
    switch (e->binop.op) {
        case OP_ADD: return l + r;
        case OP_SUB: return l - r;
        case OP_MUL: return l * r;
        case OP_DIV: if (r == 0.0) *ok = false; return l / r;
        case OP_POW: return pow(l, r);
    }
    return 0.0;
}

// S(alpha, x) for complex alpha
static double complex eval_tx(const Expr* e, double complex alpha, double x) {
    double complex sum = 0, xp = 1;
    for (int k = 0; k <= e->rootsum.s_deg; k++, xp *= x) {
        double complex c = 0, ap = 1;
        int deg = polynomial_degree(e->rootsum.s[k]);
        for (int l = 0; l <= deg; l++, ap *= alpha) {
            Rational r;
            rat_init(&r);
            poly_coeff(&r, e->rootsum.s[k], l);
            c += rat_to_double(&r) * ap;
            rat_clear(&r);
        }
        sum += c * xp;
    }
    return sum;
}

// real part of Σ alpha·log S(alpha, x) over the (numerically found) roots of q
static double eval_rootsum(const Expr* e, double x, bool* ok) {
    int n = polynomial_degree(e->rootsum.q);
    double complex* roots = malloc((size_t)n * sizeof(double complex));
    approx_complex_roots(e->rootsum.q, roots);
    double complex total = 0;
    for (int i = 0; i < n; i++) {
        double complex s = eval_tx(e, roots[i], x);
        if (s == 0) *ok = false;
        total += roots[i] * clog(s);
    }
    free(roots);
    return creal(total);
}

double expr_eval(const Expr* e, double x, bool* ok) {
    switch (e->type) {
        case EXPR_CONST: return rat_to_double(&e->constant);
        case EXPR_VAR: return x;
        case EXPR_BINOP: return eval_binop(e, x, ok);
        case EXPR_NEG: return -expr_eval(e->operand, x, ok);
        case EXPR_POLY: return poly_evaluate(e->poly, x);
        case EXPR_RATIONAL_FN: {
            double d = poly_evaluate(e->ratfn.den, x);
            if (d == 0.0) *ok = false;
            return poly_evaluate(e->ratfn.num, x) / d;
        }
        case EXPR_LN: {
            double v = fabs(expr_eval(e->ln_arg, x, ok));
            if (v == 0.0) *ok = false;
            return log(v);
        }
        case EXPR_ATAN: return atan(expr_eval(e->fn_arg, x, ok));
        case EXPR_EXP: return exp(expr_eval(e->fn_arg, x, ok));
        case EXPR_SIN: return sin(expr_eval(e->fn_arg, x, ok));
        case EXPR_COS: return cos(expr_eval(e->fn_arg, x, ok));
        case EXPR_SQRT: {
            double v = expr_eval(e->fn_arg, x, ok);
            if (v < 0.0) *ok = false;
            return sqrt(v);
        }
        case EXPR_ROOTSUM: return eval_rootsum(e, x, ok);
        case EXPR_TRIG: {
            double a = rat_to_double(&e->trig.angle) * M_PI;
            switch (e->trig.kind) {
                case TRIG_COS: return cos(a);
                case TRIG_SIN: return sin(a);
                case TRIG_COT: return cos(a) / sin(a);
                case TRIG_CSC: return 1.0 / sin(a);
                case TRIG_NONE: return a;
            }
            break;
        }
    }
    *ok = false;
    return 0.0;
}

// acc + e with e's own sum spliced in term by term (so each term's sign prints as a - b);
// takes both; acc may be NULL
Expr* expr_append_sum(Expr* acc, Expr* e) {
    if (!acc) return e;
    if (e->type == EXPR_BINOP && (e->binop.op == OP_ADD || e->binop.op == OP_SUB)) {
        Expr *l = e->binop.left, *r = e->binop.right;
        BinOp op = e->binop.op;
        e->binop.left = e->binop.right = NULL;
        free_expr(e);
        acc = expr_append_sum(acc, l);
        return op == OP_ADD ? expr_append_sum(acc, r) : create_sub(acc, r);
    }
    return create_add(acc, e);
}
