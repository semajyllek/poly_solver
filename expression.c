#include "expression.h"
#include "polynomial.h"
#include "factor.h"
#include "partial_frac.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// --- constructors ---

Expr* create_const_rat(const Rational* value) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_CONST;
    rat_init_set(&e->constant, value);
    return e;
}

Expr* create_const(double value) {
    long int_val = (long)value;
    if ((double)int_val == value) {
        Rational r = rat_from_int(int_val);
        Expr* e = create_const_rat(&r);
        rat_clear(&r);
        return e;
    }
    Rational r = rat_from_si((long)(value * 1000000), 1000000);
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
        case EXPR_CONST:
            rat_clear(&expr->constant);
            break;
    }
    free(expr);
}

// --- output ---

typedef struct { char* buf; size_t len; size_t cap; } StrBuf;

static StrBuf strbuf_new(void) {
    StrBuf sb = { malloc(64), 0, 64 };
    sb.buf[0] = '\0';
    return sb;
}

static void strbuf_append(StrBuf* sb, const char* s) {
    size_t slen = strlen(s);
    while (sb->len + slen + 1 > sb->cap) {
        sb->cap *= 2;
        sb->buf = realloc(sb->buf, sb->cap);
    }
    memcpy(sb->buf + sb->len, s, slen + 1);
    sb->len += slen;
}

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

static void expr_to_string_impl(const Expr* expr, StrBuf* sb);

// emit a term with its sign flipped: -1/2 * ln|x| → 1/2 * ln|x|
static void emit_negated_term(const Expr* expr, StrBuf* sb) {
    if (is_neg_const(expr)) {
        Rational pos; rat_init(&pos);
        rat_neg(&pos, &expr->constant);
        char* s = rat_to_string(&pos);
        strbuf_append(sb, s);
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
            strbuf_append(sb, s);
            free(s);
            strbuf_append(sb, " * ");
            expr_to_string_impl(expr->binop.right, sb);
        }
        rat_clear(&pos);
        return;
    }
    expr_to_string_impl(expr, sb);
}

// emit a + b, absorbing negative right-hand sides as subtraction
static void emit_add(const Expr* left, const Expr* right, StrBuf* sb) {
    expr_to_string_impl(left, sb);
    if (is_negative_term(right)) {
        strbuf_append(sb, " - ");
        emit_negated_term(right, sb);
    } else {
        strbuf_append(sb, " + ");
        expr_to_string_impl(right, sb);
    }
}

// emit a * b, eliding coefficients of 1 and -1
static void emit_mul(const Expr* left, const Expr* right, StrBuf* sb) {
    if (is_one_const(left)) {
        expr_to_string_impl(right, sb);
    } else if (is_neg_one_const(left)) {
        strbuf_append(sb, "-");
        expr_to_string_impl(right, sb);
    } else {
        expr_to_string_impl(left, sb);
        strbuf_append(sb, " * ");
        expr_to_string_impl(right, sb);
    }
}

// emit a / b with minimal parenthesization
static void emit_div(const Expr* left, const Expr* right, StrBuf* sb) {
    bool lp = (left->type == EXPR_BINOP &&
               (left->binop.op == OP_ADD || left->binop.op == OP_SUB));
    bool rp = (right->type == EXPR_BINOP);
    if (lp) strbuf_append(sb, "(");
    expr_to_string_impl(left, sb);
    if (lp) strbuf_append(sb, ")");
    strbuf_append(sb, " / ");
    if (rp) strbuf_append(sb, "(");
    expr_to_string_impl(right, sb);
    if (rp) strbuf_append(sb, ")");
}

// emit base^exp, rendering negative exponents as fractions
static void emit_pow(const Expr* base, const Expr* exp, StrBuf* sb) {
    if (is_neg_const(exp)) {
        Rational pos_exp; rat_init(&pos_exp);
        rat_neg(&pos_exp, &exp->constant);
        strbuf_append(sb, "1/(");
        expr_to_string_impl(base, sb);
        strbuf_append(sb, ")");
        if (!rat_is_one(&pos_exp)) {
            strbuf_append(sb, "^");
            char* s = rat_to_string(&pos_exp);
            strbuf_append(sb, s);
            free(s);
        }
        rat_clear(&pos_exp);
        return;
    }
    bool needs_parens = (base->type == EXPR_BINOP || base->type == EXPR_POLY);
    if (needs_parens) strbuf_append(sb, "(");
    expr_to_string_impl(base, sb);
    if (needs_parens) strbuf_append(sb, ")");
    strbuf_append(sb, "^");
    expr_to_string_impl(exp, sb);
}

static void expr_to_string_impl(const Expr* expr, StrBuf* sb) {
    if (!expr) { strbuf_append(sb, "NULL"); return; }

    switch (expr->type) {
        case EXPR_CONST: {
            char* s = rat_to_string(&expr->constant);
            strbuf_append(sb, s);
            free(s);
            break;
        }
        case EXPR_VAR:
            strbuf_append(sb, expr->var_name);
            break;
        case EXPR_NEG:
            strbuf_append(sb, "-(");
            expr_to_string_impl(expr->operand, sb);
            strbuf_append(sb, ")");
            break;
        case EXPR_BINOP:
            switch (expr->binop.op) {
                case OP_ADD: emit_add(expr->binop.left, expr->binop.right, sb); break;
                case OP_SUB:
                    expr_to_string_impl(expr->binop.left, sb);
                    strbuf_append(sb, " - ");
                    expr_to_string_impl(expr->binop.right, sb);
                    break;
                case OP_MUL: emit_mul(expr->binop.left, expr->binop.right, sb); break;
                case OP_DIV: emit_div(expr->binop.left, expr->binop.right, sb); break;
                case OP_POW: emit_pow(expr->binop.left, expr->binop.right, sb); break;
            }
            break;
        case EXPR_POLY: {
            char* s = poly_to_string(expr->poly);
            strbuf_append(sb, s);
            free(s);
            break;
        }
        case EXPR_RATIONAL_FN: {
            strbuf_append(sb, "(");
            char* ns = poly_to_string(expr->ratfn.num);
            strbuf_append(sb, ns);
            free(ns);
            strbuf_append(sb, ") / (");
            char* ds = poly_to_string(expr->ratfn.den);
            strbuf_append(sb, ds);
            free(ds);
            strbuf_append(sb, ")");
            break;
        }
        case EXPR_LN:
            strbuf_append(sb, "ln|");
            expr_to_string_impl(expr->ln_arg, sb);
            strbuf_append(sb, "|");
            break;
    }
}

char* expr_to_string(const Expr* expr) {
    StrBuf sb = strbuf_new();
    expr_to_string_impl(expr, &sb);
    return sb.buf;
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
    if (den->head != NULL && den->head->next == NULL &&
        den->head->exponent == 0 && rat_is_one(&den->head->coeff)) {
        free_polynomial(den);
        return wrap_poly(num);
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
        case OP_POW: {
            if (rp->head && rp->head->exponent == 0 && rp->head->next == NULL) {
                int exp = (int)rat_to_double(&rp->head->coeff);
                free_polynomial(rp);
                if (exp == 0) {
                    free_polynomial(lp);
                    Polynomial* one = create_polynomial();
                    Rational one_r = rat_one_val();
                    add_term(one, &one_r, 0);
                    rat_clear(&one_r);
                    return wrap_poly(one);
                }
                if (exp < 0) {
                    Polynomial* num = create_polynomial();
                    Rational one_r = rat_one_val();
                    add_term(num, &one_r, 0);
                    rat_clear(&one_r);
                    Polynomial* den = poly_copy(lp);
                    for (int i = 1; i < -exp; i++) {
                        Polynomial* tmp = poly_multiply(den, lp);
                        free_polynomial(den);
                        den = tmp;
                    }
                    free_polynomial(lp);
                    return wrap_ratfn(num, den);
                }
                Polynomial* result = poly_copy(lp);
                for (int i = 1; i < exp; i++) {
                    Polynomial* tmp = poly_multiply(result, lp);
                    free_polynomial(result);
                    result = tmp;
                }
                free_polynomial(lp);
                return wrap_poly(result);
            }
            free_polynomial(lp); free_polynomial(rp);
            return NULL;
        }
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

static Expr* make_linear_expr(const Rational* root) {
    Polynomial* p = create_polynomial();
    Rational one_r = rat_one_val();
    add_term(p, &one_r, 1);
    rat_clear(&one_r);
    Rational neg_root;
    rat_init(&neg_root);
    rat_neg(&neg_root, root);
    add_term(p, &neg_root, 0);
    rat_clear(&neg_root);
    return create_poly_expr(p);
}

static Expr* integrate_pf_term(PFTerm* t) {
    if (t->is_quadratic) return NULL;

    // linear factor (x - r): constant term negated gives root
    Rational root;
    rat_init(&root);
    rat_neg(&root, &t->factor->head->next->coeff);

    Expr* result = NULL;

    for (int k = 0; k < t->power; k++) {
        if (rat_is_zero(&t->coeffs[k])) continue;

        int power = k + 1;
        Expr* term;

        if (power == 1) {
            term = create_mul(create_const_rat(&t->coeffs[k]),
                              create_ln(make_linear_expr(&root)));
        } else {
            Rational new_exp = rat_from_int(1 - power);
            Rational divisor = rat_from_int(1 - power);
            Rational new_coeff;
            rat_init(&new_coeff);
            rat_div(&new_coeff, &t->coeffs[k], &divisor);
            term = create_mul(
                create_const_rat(&new_coeff),
                create_pow(make_linear_expr(&root), create_const_rat(&new_exp)));
            rat_clear(&new_coeff);
            rat_clear(&new_exp);
            rat_clear(&divisor);
        }

        if (result == NULL) result = term;
        else result = create_add(result, term);
    }

    rat_clear(&root);
    return result;
}

Expr* integrate(Expr* expr, const char* variable) {
    if (!expr) return NULL;
    (void)variable;

    Expr* canonical = canonicalize(expr);
    if (!canonical) return NULL;

    if (canonical->type == EXPR_POLY) {
        Polynomial* integral = poly_integral(canonical->poly);
        free_expr(canonical);
        return wrap_poly(integral);
    }

    if (canonical->type == EXPR_RATIONAL_FN) {
        PartialFractionDecomp* pf = partial_fractions(
            canonical->ratfn.num, canonical->ratfn.den);
        free_expr(canonical);
        if (!pf) return NULL;

        Expr* result = NULL;

        if (pf->poly_part && !poly_is_zero(pf->poly_part)) {
            Polynomial* poly_int = poly_integral(pf->poly_part);
            result = wrap_poly(poly_int);
        }

        for (int i = 0; i < pf->num_terms; i++) {
            Expr* term_integral = integrate_pf_term(&pf->terms[i]);
            if (!term_integral) continue;
            if (result == NULL) result = term_integral;
            else result = create_add(result, term_integral);
        }

        free_pf_decomp(pf);
        return result;
    }

    free_expr(canonical);
    return NULL;
}
