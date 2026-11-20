#include "polynomial.h"
#include "strbuf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Term* create_term(const Rational* coeff, int exponent) {
    Term* term = malloc(sizeof(Term));
    rat_init_set(&term->coeff, coeff);
    term->exponent = exponent;
    term->var_id = 0;
    term->next = NULL;
    return term;
}

Polynomial* create_polynomial(void) {
    Polynomial* poly = malloc(sizeof(Polynomial));
    poly->head = NULL;
    return poly;
}

// builds a polynomial term by term in descending exponent order, in O(1) per term
typedef struct {
    Polynomial* p;
    Term** tail;
} Appender;

static Appender app_new(void) {
    Polynomial* p = create_polynomial();
    return (Appender){ p, &p->head };
}

// appends c·x^e (e below every exponent so far); zeros are skipped
static void app_push(Appender* a, const mpq_t c, int e) {
    if (mpq_sgn(c) == 0) return;
    Term* t = malloc(sizeof(Term));
    mpq_init(t->coeff.val);
    mpq_set(t->coeff.val, c);
    t->exponent = e;
    t->var_id = 0;
    t->next = NULL;
    *a->tail = t;
    a->tail = &t->next;
}

// lowest exponent of a nonzero polynomial (its last term)
static int low_exponent(const Polynomial* p) {
    const Term* t = p->head;
    while (t->next) t = t->next;
    return t->exponent;
}

// dense coefficient window c[0..n), c[i] the coefficient of x^(lo + i)
typedef struct {
    mpq_t* c;
    int lo, n;
} Dense;

static Dense dense_zero(int lo, int hi) {
    Dense d = { malloc(sizeof(mpq_t) * (size_t)(hi - lo + 1)), lo, hi - lo + 1 };
    for (int i = 0; i < d.n; i++) mpq_init(d.c[i]);
    return d;
}

static void dense_add_poly(Dense* d, const Polynomial* p) {
    for (const Term* t = p->head; t; t = t->next) mpq_add(d->c[t->exponent - d->lo], d->c[t->exponent - d->lo], t->coeff.val);
}

// the nonzero entries as a polynomial; frees the window
static Polynomial* dense_take(Dense* d) {
    Appender a = app_new();
    for (int i = d->n - 1; i >= 0; i--) { app_push(&a, d->c[i], d->lo + i); mpq_clear(d->c[i]); }
    free(d->c);
    return a.p;
}

void add_term(Polynomial* poly, const Rational* coeff, int exponent) {
    if (rat_is_zero(coeff)) return;

    Term* new_term = create_term(coeff, exponent);
    if (poly->head == NULL) {
        poly->head = new_term;
        return;
    }

    Term* current = poly->head;
    Term* previous = NULL;

    while (current != NULL && current->exponent > exponent) {
        previous = current;
        current = current->next;
    }

    if (current != NULL && current->exponent == exponent) {
        rat_add(&current->coeff, &current->coeff, coeff);
        rat_clear(&new_term->coeff);
        free(new_term);
        if (rat_is_zero(&current->coeff)) {
            if (previous == NULL)
                poly->head = current->next;
            else
                previous->next = current->next;
            rat_clear(&current->coeff);
            free(current);
        }
    } else {
        new_term->next = current;
        if (previous == NULL)
            poly->head = new_term;
        else
            previous->next = new_term;
    }
}

Polynomial* poly_copy(const Polynomial* poly) {
    Appender a = app_new();
    for (const Term* t = poly->head; t; t = t->next) app_push(&a, t->coeff.val, t->exponent);
    return a.p;
}

bool poly_is_zero(const Polynomial* poly) {
    return poly->head == NULL;
}

int polynomial_degree(const Polynomial* poly) {
    if (poly->head == NULL) return -1;
    return poly->head->exponent;
}

// p ± q by merging the two sorted term lists
static Polynomial* merge(const Polynomial* p, const Polynomial* q, int sign) {
    Appender a = app_new();
    const Term *s = p->head, *t = q->head;
    mpq_t c;
    mpq_init(c);
    while (s || t) {
        if (t && (!s || t->exponent > s->exponent)) {            // only q has this exponent
            if (sign > 0) mpq_set(c, t->coeff.val); else mpq_neg(c, t->coeff.val);
            app_push(&a, c, t->exponent);
            t = t->next;
        } else if (!t || s->exponent > t->exponent) {            // only p has it
            app_push(&a, s->coeff.val, s->exponent);
            s = s->next;
        } else {                                                 // both: combine
            if (sign > 0) mpq_add(c, s->coeff.val, t->coeff.val); else mpq_sub(c, s->coeff.val, t->coeff.val);
            app_push(&a, c, s->exponent);
            s = s->next; t = t->next;
        }
    }
    mpq_clear(c);
    return a.p;
}

Polynomial* poly_add(const Polynomial* poly1, const Polynomial* poly2) { return merge(poly1, poly2, 1); }

Polynomial* poly_subtract(const Polynomial* poly1, const Polynomial* poly2) { return merge(poly1, poly2, -1); }

static int term_count(const Polynomial* p) {
    int n = 0;
    for (const Term* t = p->head; t; t = t->next) n++;
    return n;
}

// very sparse factors (x^1000·y): insert each product term instead of a huge window
static Polynomial* sparse_multiply(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* result = create_polynomial();
    Rational prod;
    rat_init(&prod);
    for (const Term* s = poly1->head; s; s = s->next)
        for (const Term* t = poly2->head; t; t = t->next) {
            rat_mul(&prod, &s->coeff, &t->coeff);
            add_term(result, &prod, s->exponent + t->exponent);
        }
    rat_clear(&prod);
    return result;
}

// schoolbook product accumulated in a dense window over [lo1 + lo2, deg1 + deg2]
Polynomial* poly_multiply(const Polynomial* poly1, const Polynomial* poly2) {
    if (!poly1->head || !poly2->head) return create_polynomial();
    long span = (long)poly1->head->exponent + poly2->head->exponent - low_exponent(poly1) - low_exponent(poly2);
    if (span > 64L * term_count(poly1) * term_count(poly2) + 1024) return sparse_multiply(poly1, poly2);
    Dense d = dense_zero(low_exponent(poly1) + low_exponent(poly2), poly1->head->exponent + poly2->head->exponent);
    mpq_t prod;
    mpq_init(prod);
    for (const Term* s = poly1->head; s; s = s->next)
        for (const Term* t = poly2->head; t; t = t->next) {
            mpq_mul(prod, s->coeff.val, t->coeff.val);
            mpq_add(d.c[s->exponent + t->exponent - d.lo], d.c[s->exponent + t->exponent - d.lo], prod);
        }
    mpq_clear(prod);
    return dense_take(&d);
}

// division in place on a dense copy of the dividend: each step clears the top coefficient
PolyDivResult poly_divmod(const Polynomial* dividend, const Polynomial* divisor) {
    if (divisor == NULL || divisor->head == NULL) {
        fprintf(stderr, "Error: division by zero polynomial.\n");
        return (PolyDivResult){ create_polynomial(), create_polynomial() };
    }
    int db = divisor->head->exponent, da = polynomial_degree(dividend);
    if (da < db) return (PolyDivResult){ create_polynomial(), poly_copy(dividend) };
    int lo = low_exponent(dividend), lb = low_exponent(divisor);
    Dense r = dense_zero(lo < lb ? lo : lb, da), q = dense_zero(0, da - db);
    dense_add_poly(&r, dividend);
    mpq_t inv, prod;
    mpq_inits(inv, prod, NULL);
    mpq_inv(inv, divisor->head->coeff.val);
    for (int e = da; e >= db; e--) {
        mpq_t* top = &r.c[e - r.lo];
        if (mpq_sgn(*top) == 0) continue;
        mpq_t* qc = &q.c[e - db - q.lo];
        mpq_mul(*qc, *top, inv);                                 // next quotient term
        for (const Term* t = divisor->head->next; t; t = t->next) {
            mpq_mul(prod, *qc, t->coeff.val);
            mpq_sub(r.c[e - db + t->exponent - r.lo], r.c[e - db + t->exponent - r.lo], prod);
        }
        mpq_set_ui(*top, 0, 1);
    }
    mpq_clears(inv, prod, NULL);
    return (PolyDivResult){ dense_take(&q), dense_take(&r) };
}

Polynomial* poly_divide(const Polynomial* dividend, const Polynomial* divisor) {
    PolyDivResult dr = poly_divmod(dividend, divisor);
    free_polynomial(dr.remainder);
    return dr.quotient;
}

Polynomial* poly_mod(const Polynomial* dividend, const Polynomial* divisor) {
    PolyDivResult dr = poly_divmod(dividend, divisor);
    free_polynomial(dr.quotient);
    return dr.remainder;
}

void poly_make_monic(Polynomial* poly) {
    if (poly->head == NULL) return;
    if (rat_is_one(&poly->head->coeff)) return;
    mpq_t inv;
    mpq_init(inv);
    mpq_inv(inv, poly->head->coeff.val);
    for (Term* t = poly->head; t; t = t->next) mpq_mul(t->coeff.val, t->coeff.val, inv);
    mpq_clear(inv);
}

Polynomial* poly_gcd(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* r0 = poly_copy(poly1);
    Polynomial* r1 = poly_copy(poly2);

    poly_make_monic(r1);
    while (!poly_is_zero(r1)) {
        Polynomial* r2 = poly_mod(r0, r1);
        poly_make_monic(r2);                                     // monic remainders: smaller coefficients
        free_polynomial(r0);
        r0 = r1;
        r1 = r2;
    }

    free_polynomial(r1);
    poly_make_monic(r0);
    return r0;
}

Polynomial* poly_lcm(const Polynomial* poly1, const Polynomial* poly2) {
    Polynomial* product = poly_multiply(poly1, poly2);
    Polynomial* gcd = poly_gcd(poly1, poly2);
    Polynomial* lcm = poly_divide(product, gcd);
    free_polynomial(product);
    free_polynomial(gcd);
    return lcm;
}

Polynomial* poly_constant(const Rational* c) {
    Polynomial* p = create_polynomial();
    add_term(p, c, 0);
    return p;
}

Polynomial* poly_pow(const Polynomial* base, unsigned exp) {
    Rational one = rat_one_val();
    Polynomial* result = poly_constant(&one);
    rat_clear(&one);
    Polynomial* sq = poly_copy(base);
    while (exp > 0) {
        if (exp & 1u) { // fold the current square into the result
            Polynomial* t = poly_multiply(result, sq);
            free_polynomial(result);
            result = t;
        }
        exp >>= 1;
        if (exp) {
            Polynomial* t = poly_multiply(sq, sq);
            free_polynomial(sq);
            sq = t;
        }
    }
    free_polynomial(sq);
    return result;
}

Polynomial* poly_scale(const Polynomial* p, const Rational* c) {
    Appender a = app_new();
    mpq_t t;
    mpq_init(t);
    for (const Term* term = p->head; term; term = term->next) {
        mpq_mul(t, term->coeff.val, c->val);
        app_push(&a, t, term->exponent);
    }
    mpq_clear(t);
    return a.p;
}

// next = old - q * cur (one row of the extended Euclidean recurrence)
static Polynomial* euclid_row(const Polynomial* old, const Polynomial* q, const Polynomial* cur) {
    Polynomial* qc = poly_multiply(q, cur);
    Polynomial* r = poly_subtract(old, qc);
    free_polynomial(qc);
    return r;
}

PolyBezout poly_xgcd(const Polynomial* a, const Polynomial* b) {
    Rational one = rat_one_val();
    Polynomial *r0 = poly_copy(a), *r1 = poly_copy(b);
    Polynomial *s0 = poly_constant(&one), *s1 = create_polynomial();
    Polynomial *t0 = create_polynomial(), *t1 = poly_constant(&one);
    while (!poly_is_zero(r1)) {
        PolyDivResult dr = poly_divmod(r0, r1);
        Polynomial* s2 = euclid_row(s0, dr.quotient, s1);
        Polynomial* t2 = euclid_row(t0, dr.quotient, t1);
        free_polynomial(dr.quotient); free_polynomial(r0); free_polynomial(s0); free_polynomial(t0);
        r0 = r1; r1 = dr.remainder;
        s0 = s1; s1 = s2;
        t0 = t1; t1 = t2;
    }
    // normalize so the gcd is monic
    Rational inv;
    rat_init(&inv);
    if (r0->head) mpq_inv(inv.val, r0->head->coeff.val);
    else mpq_set_si(inv.val, 1, 1);
    PolyBezout bz = { poly_scale(s0, &inv), poly_scale(t0, &inv), poly_scale(r0, &inv) };
    free_polynomial(r0); free_polynomial(r1); free_polynomial(s0);
    free_polynomial(s1); free_polynomial(t0); free_polynomial(t1);
    rat_clear(&inv); rat_clear(&one);
    return bz;
}

void free_bezout(PolyBezout* bz) {
    free_polynomial(bz->s);
    free_polynomial(bz->t);
    free_polynomial(bz->g);
}

void poly_coeff(Rational* out, const Polynomial* poly, int degree) {
    for (Term* t = poly->head; t; t = t->next) {
        if (t->exponent == degree) { rat_set(out, &t->coeff); return; }
        if (t->exponent < degree) break; // terms are sorted descending
    }
    mpq_set_si(out->val, 0, 1);
}

void poly_eval_rat(Rational* out, const Polynomial* poly, const Rational* x) {
    mpq_set_si(out->val, 0, 1);
    int prev = poly->head ? poly->head->exponent : 0;
    for (Term* t = poly->head; t; t = t->next) {
        for (int i = t->exponent; i < prev; i++) rat_mul(out, out, x); // skip the gap
        rat_add(out, out, &t->coeff);
        prev = t->exponent;
    }
    for (int i = 0; i < prev; i++) rat_mul(out, out, x);
}

DPoly dpoly_from(const Polynomial* p) {
    DPoly d = { NULL, polynomial_degree(p) };
    if (d.deg < 0) return d;
    d.c = calloc((size_t)d.deg + 1, sizeof(double));
    for (const Term* t = p->head; t; t = t->next) d.c[t->exponent] = mpq_get_d(t->coeff.val);
    return d;
}

// Horner; a zero coefficient adds 0 exactly, so this matches poly_evaluate bit for bit
double dpoly_eval(const DPoly* p, double x) {
    if (p->deg < 0) return 0.0;
    double r = p->c[p->deg];
    for (int i = p->deg - 1; i >= 0; i--) r = r * x + p->c[i];
    return r;
}

void dpoly_free(DPoly* p) { free(p->c); p->c = NULL; }

double poly_evaluate(const Polynomial* poly, double x) {
    if (poly->head == NULL) return 0.0;

    Term* current = poly->head;
    double result = rat_to_double(&current->coeff);
    int prev_exp = current->exponent;
    current = current->next;

    while (current != NULL) {
        int gap = prev_exp - current->exponent;
        for (int i = 0; i < gap; i++)
            result *= x;
        result += rat_to_double(&current->coeff);
        prev_exp = current->exponent;
        current = current->next;
    }

    for (int i = 0; i < prev_exp; i++)
        result *= x;

    return result;
}

Polynomial* poly_derivative(const Polynomial* poly) {
    Appender a = app_new();
    mpq_t c;
    mpq_init(c);
    for (const Term* t = poly->head; t && t->exponent > 0; t = t->next) {
        mpq_set_si(c, t->exponent, 1);
        mpq_mul(c, c, t->coeff.val);
        app_push(&a, c, t->exponent - 1);
    }
    mpq_clear(c);
    return a.p;
}

Polynomial* poly_integral(const Polynomial* poly) {
    Appender a = app_new();
    mpq_t c;
    mpq_init(c);
    for (const Term* t = poly->head; t; t = t->next) {
        mpq_set_si(c, t->exponent + 1, 1);
        mpq_div(c, t->coeff.val, c);
        app_push(&a, c, t->exponent + 1);
    }
    mpq_clear(c);
    return a.p;
}

double poly_definite_integral(const Polynomial* poly, double lower, double upper) {
    Polynomial* integral = poly_integral(poly);
    double result = poly_evaluate(integral, upper) - poly_evaluate(integral, lower);
    free_polynomial(integral);
    return result;
}

void free_polynomial(Polynomial* poly) {
    if (!poly) return;
    Term* current = poly->head;
    while (current != NULL) {
        Term* temp = current;
        current = current->next;
        rat_clear(&temp->coeff);
        free(temp);
    }
    free(poly);
}

void print_polynomial(const Polynomial* poly) {
    char* s = poly_to_string(poly);
    printf("%s", s);
    free(s);
}

typedef struct {
    const Rational* abs_c; // |coefficient|
    int e;                 // exponent
    const char* var;       // variable name
} TermFmt;

// append one term (sign handled by caller) as "c", "x^e", or "cx^e"
static void append_term_body(StrBuf* sb, const TermFmt* t) {
    if (t->e == 0) { sb_append_owned(sb, rat_to_string(t->abs_c)); return; }
    if (!rat_is_one(t->abs_c)) sb_append_owned(sb, rat_to_string(t->abs_c));
    sb_append(sb, t->var);
    if (t->e != 1) sb_appendf(sb, "^%d", t->e);
}

char* poly_to_string_var(const Polynomial* poly, const char* var) {
    StrBuf sb = sb_new();
    if (poly->head == NULL) { sb_append(&sb, "0"); return sb_finish(&sb); }
    Rational c;
    rat_init(&c);
    for (Term* t = poly->head; t; t = t->next) {
        bool neg = rat_is_negative(&t->coeff);
        if (t != poly->head) sb_append(&sb, neg ? " - " : " + ");
        else if (neg) sb_append(&sb, "-");
        rat_abs(&c, &t->coeff);
        TermFmt fmt = { &c, t->exponent, var };
        append_term_body(&sb, &fmt);
    }
    rat_clear(&c);
    return sb_finish(&sb);
}

char* poly_to_string(const Polynomial* poly) {
    return poly_to_string_var(poly, "x");
}

QuadRoots poly_quad_roots(const Polynomial* f) {
    QuadRoots q;
    Rational a, b, c, d, t;
    rat_init(&q.p); rat_init(&q.r); mpz_init(q.m);
    rat_init(&a); rat_init(&b); rat_init(&c); rat_init(&d); rat_init(&t);
    poly_coeff(&a, f, 2); poly_coeff(&b, f, 1); poly_coeff(&c, f, 0);
    rat_mul(&d, &b, &b);                       // D = b^2 - 4ac
    rat_mul(&t, &a, &c);
    mpq_mul_2exp(t.val, t.val, 2);
    rat_sub(&d, &d, &t);
    q.is_complex = rat_is_negative(&d);
    rat_abs(&d, &d);
    rat_mul(&t, &a, &a);                       // |D| / (4a^2)
    mpq_mul_2exp(t.val, t.val, 2);
    rat_div(&d, &d, &t);
    rat_sqrt_split(&d, &q.r, q.m);
    rat_div(&q.p, &b, &a);                     // p = -b / 2a
    mpq_div_2exp(q.p.val, q.p.val, 1);
    rat_neg(&q.p, &q.p);
    rat_clear(&a); rat_clear(&b); rat_clear(&c); rat_clear(&d); rat_clear(&t);
    return q;
}

void quad_roots_clear(QuadRoots* q) {
    rat_clear(&q->p); rat_clear(&q->r); mpz_clear(q->m);
}
