#include "expression.h"
#include "polynomial.h"
#include "factor.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// --- constructors ---

Expr* create_const_rat(Rational value) {
    Expr* e = malloc(sizeof(Expr));
    e->type = EXPR_CONST;
    e->constant = value;
    return e;
}

Expr* create_const(double value) {
    // convert double to rational: handle integers exactly, fractions approximately
    long long int_val = (long long)value;
    if ((double)int_val == value)
        return create_const_rat(rat_from_int(int_val));
    // fallback: use large denominator for non-integer doubles
    return create_const_rat(rat_create((long long)(value * 1000000), 1000000));
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
        case EXPR_CONST:
            break;
    }
    free(expr);
}

// --- output ---

// helper: append to a dynamically-growing string buffer
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

static void expr_to_string_impl(const Expr* expr, StrBuf* sb) {
    if (!expr) { strbuf_append(sb, "NULL"); return; }

    switch (expr->type) {
        case EXPR_CONST: {
            char* s = rat_to_string(expr->constant);
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
        case EXPR_BINOP: {
            bool needs_parens = (expr->binop.op == OP_ADD || expr->binop.op == OP_SUB);
            if (needs_parens) strbuf_append(sb, "(");
            expr_to_string_impl(expr->binop.left, sb);
            switch (expr->binop.op) {
                case OP_ADD: strbuf_append(sb, " + "); break;
                case OP_SUB: strbuf_append(sb, " - "); break;
                case OP_MUL: strbuf_append(sb, " * "); break;
                case OP_DIV: strbuf_append(sb, " / "); break;
                case OP_POW: strbuf_append(sb, "^"); break;
            }
            expr_to_string_impl(expr->binop.right, sb);
            if (needs_parens) strbuf_append(sb, ")");
            break;
        }
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
// converts an expression tree into EXPR_POLY or EXPR_RATIONAL_FN

// helper: wrap a polynomial in an Expr (takes ownership of poly)
static Expr* wrap_poly(Polynomial* poly) {
    return create_poly_expr(poly);
}

// helper: wrap a rational function, simplifying if denominator is 1
static Expr* wrap_ratfn(Polynomial* num, Polynomial* den) {
    // if denominator is a nonzero constant with value 1, return just the numerator
    if (den->head != NULL && den->head->next == NULL &&
        den->head->exponent == 0 && rat_is_one(den->head->coeff)) {
        free_polynomial(den);
        return wrap_poly(num);
    }
    // reduce by GCD
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
    // if denominator is now constant 1
    if (den->head != NULL && den->head->next == NULL &&
        den->head->exponent == 0 && rat_is_one(den->head->coeff)) {
        free_polynomial(den);
        return wrap_poly(num);
    }
    return create_rational_fn(num, den);
}

// canonicalize a polynomial expr with a polynomial expr via binary op
static Expr* canon_poly_op(BinOp op, Polynomial* lp, Polynomial* rp) {
    switch (op) {
        case OP_ADD: {
            Polynomial* result = poly_add(lp, rp);
            free_polynomial(lp);
            free_polynomial(rp);
            return wrap_poly(result);
        }
        case OP_SUB: {
            Polynomial* result = poly_subtract(lp, rp);
            free_polynomial(lp);
            free_polynomial(rp);
            return wrap_poly(result);
        }
        case OP_MUL: {
            Polynomial* result = poly_multiply(lp, rp);
            free_polynomial(lp);
            free_polynomial(rp);
            return wrap_poly(result);
        }
        case OP_DIV:
            return wrap_ratfn(lp, rp);
        case OP_POW: {
            // only support integer exponent
            if (rp->head && rp->head->exponent == 0 && rp->head->next == NULL) {
                int exp = (int)rat_to_double(rp->head->coeff);
                free_polynomial(rp);
                if (exp == 0) {
                    free_polynomial(lp);
                    Polynomial* one = create_polynomial();
                    add_term(one, rat_one(), 0);
                    return wrap_poly(one);
                }
                if (exp < 0) {
                    // 1 / poly^|exp|
                    Polynomial* num = create_polynomial();
                    add_term(num, rat_one(), 0);
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
            // non-constant exponent: can't canonicalize
            free_polynomial(lp);
            free_polynomial(rp);
            return NULL;
        }
    }
    free_polynomial(lp);
    free_polynomial(rp);
    return NULL;
}

// extract polynomial from canonical expr, caller must free the returned poly
static Polynomial* extract_poly(Expr* e) {
    if (e->type == EXPR_POLY) {
        Polynomial* p = e->poly;
        e->poly = NULL; // transfer ownership
        return p;
    }
    return NULL;
}

// extract num/den from canonical expr
static bool extract_ratfn(Expr* e, Polynomial** num, Polynomial** den) {
    if (e->type == EXPR_POLY) {
        *num = e->poly;
        e->poly = NULL;
        *den = create_polynomial();
        add_term(*den, rat_one(), 0);
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
            add_term(p, expr->constant, 0);
            return wrap_poly(p);
        }
        case EXPR_VAR: {
            Polynomial* p = create_polynomial();
            add_term(p, rat_one(), 1);
            return wrap_poly(p);
        }
        case EXPR_NEG: {
            Expr* inner = canonicalize(expr->operand);
            if (!inner) return NULL;
            if (inner->type == EXPR_POLY) {
                // negate all coefficients
                Polynomial* neg = create_polynomial();
                Term* t = inner->poly->head;
                while (t) {
                    add_term(neg, rat_neg(t->coeff), t->exponent);
                    t = t->next;
                }
                free_expr(inner);
                return wrap_poly(neg);
            }
            if (inner->type == EXPR_RATIONAL_FN) {
                Polynomial* neg_num = create_polynomial();
                Term* t = inner->ratfn.num->head;
                while (t) {
                    add_term(neg_num, rat_neg(t->coeff), t->exponent);
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
                free_expr(left);
                free_expr(right);
                return NULL;
            }

            // if both are polynomials, operate directly
            if (left->type == EXPR_POLY && right->type == EXPR_POLY) {
                Polynomial* lp = extract_poly(left);
                Polynomial* rp = extract_poly(right);
                free(left);
                free(right);
                return canon_poly_op(expr->binop.op, lp, rp);
            }

            // lift to rational functions
            Polynomial *ln, *ld, *rn, *rd;
            if (!extract_ratfn(left, &ln, &ld) || !extract_ratfn(right, &rn, &rd)) {
                free_expr(left);
                free_expr(right);
                return NULL;
            }
            free(left);
            free(right);

            switch (expr->binop.op) {
                case OP_ADD: case OP_SUB: {
                    // a/b +- c/d = (a*d +- c*b) / (b*d)
                    Polynomial* ad = poly_multiply(ln, rd);
                    Polynomial* cb = poly_multiply(rn, ld);
                    Polynomial* num = (expr->binop.op == OP_ADD) ?
                        poly_add(ad, cb) : poly_subtract(ad, cb);
                    Polynomial* den = poly_multiply(ld, rd);
                    free_polynomial(ad);
                    free_polynomial(cb);
                    free_polynomial(ln);
                    free_polynomial(ld);
                    free_polynomial(rn);
                    free_polynomial(rd);
                    return wrap_ratfn(num, den);
                }
                case OP_MUL: {
                    Polynomial* num = poly_multiply(ln, rn);
                    Polynomial* den = poly_multiply(ld, rd);
                    free_polynomial(ln);
                    free_polynomial(ld);
                    free_polynomial(rn);
                    free_polynomial(rd);
                    return wrap_ratfn(num, den);
                }
                case OP_DIV: {
                    Polynomial* num = poly_multiply(ln, rd);
                    Polynomial* den = poly_multiply(ld, rn);
                    free_polynomial(ln);
                    free_polynomial(ld);
                    free_polynomial(rn);
                    free_polynomial(rd);
                    return wrap_ratfn(num, den);
                }
                case OP_POW: {
                    // only if right is a constant integer
                    free_polynomial(rd);
                    Polynomial* result_num = canon_poly_op(OP_POW, ln, rn) ? NULL : NULL;
                    // simplified: handle pow of rational function
                    // for now, only support polynomial base
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
    }
    return NULL;
}

Expr* simplify(Expr* expr) {
    return canonicalize(expr);
}

// --- stubs for later phases ---

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

    // build expression: content * factor1^m1 * factor2^m2 * ...
    Expr* result = NULL;

    if (f->count == 0) {
        result = create_const_rat(f->content);
        free_factorization(f);
        return result;
    }

    for (int i = 0; i < f->count; i++) {
        Expr* factor_expr = create_poly_expr(poly_copy(f->factors[i]));
        if (f->multiplicities[i] > 1)
            factor_expr = create_pow(factor_expr, create_const_rat(rat_from_int(f->multiplicities[i])));
        if (result == NULL)
            result = factor_expr;
        else
            result = create_mul(result, factor_expr);
    }

    if (!rat_is_one(f->content))
        result = create_mul(create_const_rat(f->content), result);

    free_factorization(f);
    return result;
}

Expr* differentiate(Expr* expr, const char* variable) {
    if (!expr) return NULL;
    (void)variable; // only single-variable x for now

    Expr* canonical = canonicalize(expr);
    if (!canonical) return NULL;

    if (canonical->type == EXPR_POLY) {
        Polynomial* deriv = poly_derivative(canonical->poly);
        free_expr(canonical);
        return wrap_poly(deriv);
    }

    if (canonical->type == EXPR_RATIONAL_FN) {
        // quotient rule: (f'g - fg') / g^2
        Polynomial* f = canonical->ratfn.num;
        Polynomial* g = canonical->ratfn.den;
        Polynomial* fp = poly_derivative(f);
        Polynomial* gp = poly_derivative(g);
        Polynomial* fp_g = poly_multiply(fp, g);
        Polynomial* f_gp = poly_multiply(f, gp);
        Polynomial* num = poly_subtract(fp_g, f_gp);
        Polynomial* den = poly_multiply(g, g);
        free_polynomial(fp);
        free_polynomial(gp);
        free_polynomial(fp_g);
        free_polynomial(f_gp);
        free_expr(canonical);
        return wrap_ratfn(num, den);
    }

    free_expr(canonical);
    return NULL;
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

    // rational function integration not yet supported
    free_expr(canonical);
    return NULL;
}
