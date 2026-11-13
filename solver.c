#include "solver.h"
#include "factor.h"
#include "strbuf.h"
#include <stdlib.h>
#include <string.h>

#define DECIMAL_DIGITS 10 // digits shown for irrational roots (truncated)

// --- building the solution set ---

static void push_rational(Solutions* sol, const Rational* r) {
    sol->rational_roots = realloc(sol->rational_roots, (size_t)(sol->num_rational + 1) * sizeof(Rational));
    rat_init_set(&sol->rational_roots[sol->num_rational++], r);
}

// root -c0/c1 of a linear factor, repeated by multiplicity
static void solve_linear(Solutions* sol, const Polynomial* f, int mult) {
    Rational c0, c1, root;
    rat_init(&c0); rat_init(&c1); rat_init(&root);
    poly_coeff(&c0, f, 0);
    poly_coeff(&c1, f, 1);
    rat_div(&root, &c0, &c1);
    rat_neg(&root, &root);
    for (int m = 0; m < mult; m++) push_rational(sol, &root);
    rat_clear(&c0); rat_clear(&c1); rat_clear(&root);
}

static void midpoint(Rational* out, const RootInterval* r) {
    rat_add(out, &r->lo, &r->hi);
    mpq_div_2exp(out->val, out->val, 1);
}

// fill the legacy double / non-real-factor views for one irrational factor
static void push_legacy_views(Solutions* sol, const IrrationalFactor* e) {
    Rational mid;
    rat_init(&mid);
    for (int i = 0; i < e->num_real; i++) {
        midpoint(&mid, &e->real[i]);
        for (int m = 0; m < e->multiplicity; m++) {
            sol->irrational_roots = realloc(sol->irrational_roots, (size_t)(sol->num_irrational + 1) * sizeof(double));
            sol->irrational_roots[sol->num_irrational++] = rat_to_double(&mid);
        }
    }
    rat_clear(&mid);
    if (e->num_real == polynomial_degree(e->factor)) return; // all roots are real
    for (int m = 0; m < e->multiplicity; m++) {
        sol->irreducible_factors = realloc(sol->irreducible_factors, (size_t)(sol->num_irreducible + 1) * sizeof(Polynomial*));
        sol->irreducible_factors[sol->num_irreducible++] = poly_copy(e->factor);
    }
}

// isolate and refine the real roots of an irreducible factor of degree >= 2
static void solve_irrational(Solutions* sol, const Polynomial* f, int mult) {
    RootIntervalList L = isolate_real_roots(f);
    Rational width; // 10^-(digits+2): well below the last digit shown
    rat_init(&width);
    mpz_ui_pow_ui(mpq_denref(width.val), 10, DECIMAL_DIGITS + 2);
    mpz_set_ui(mpq_numref(width.val), 1);
    for (int i = 0; i < L.count; i++) refine_root(&L.items[i], f, &width);
    rat_clear(&width);
    sol->irr = realloc(sol->irr, (size_t)(sol->num_irr + 1) * sizeof(IrrationalFactor));
    IrrationalFactor* e = &sol->irr[sol->num_irr++];
    *e = (IrrationalFactor){ poly_copy(f), mult, L.items, L.count };
    push_legacy_views(sol, e);
}

Solutions* solve_polynomial(const Polynomial* poly) {
    Solutions* sol = calloc(1, sizeof(Solutions));
    if (poly->head == NULL) return sol;
    Factorization* f = factorize(poly);
    for (int i = 0; i < f->count; i++) {
        if (polynomial_degree(f->factors[i]) == 1) solve_linear(sol, f->factors[i], f->multiplicities[i]);
        else solve_irrational(sol, f->factors[i], f->multiplicities[i]);
    }
    free_factorization(f);
    return sol;
}

// --- output ---

// sqrt-term text "r*sqrt(m)" with the trivial parts elided ("1" if both trivial)
static void append_radical(StrBuf* sb, const Rational* r, const mpz_t m) {
    bool unit_r = rat_is_one(r), unit_m = mpz_cmp_ui(m, 1) == 0;
    if (!unit_r || unit_m) sb_append_owned(sb, rat_to_string(r));
    if (!unit_r && !unit_m) sb_append(sb, "*");
    if (!unit_m) { sb_append(sb, "sqrt("); sb_append_owned(sb, mpz_get_str(NULL, 10, m)); sb_append(sb, ")"); }
}

static void separator(StrBuf* sb, bool* first) {
    if (!*first) sb_append(sb, ", ");
    *first = false;
}

// " ~ 1.4142135623" from a refined isolating interval
static void append_approx(StrBuf* sb, const RootInterval* r) {
    Rational mid;
    rat_init(&mid);
    midpoint(&mid, r);
    sb_append(sb, " ~ ");
    sb_append_owned(sb, rat_to_decimal(&mid, DECIMAL_DIGITS));
    rat_clear(&mid);
}

static void append_multiplicity(StrBuf* sb, int mult) {
    if (mult > 1) sb_appendf(sb, " (multiplicity %d)", mult);
}

// a quadratic's roots p ± r*sqrt(m) (real) or p ± r*sqrt(m)*i (complex), exactly
typedef struct {
    Rational p, r;
    mpz_t m;
    bool complex;
} QuadRoots;

static QuadRoots quad_roots(const Polynomial* f) {
    QuadRoots q;
    Rational a, b, c, d, t;
    rat_init(&q.p); rat_init(&q.r); mpz_init(q.m);
    rat_init(&a); rat_init(&b); rat_init(&c); rat_init(&d); rat_init(&t);
    poly_coeff(&a, f, 2); poly_coeff(&b, f, 1); poly_coeff(&c, f, 0);
    rat_mul(&d, &b, &b);                       // D = b^2 - 4ac
    rat_mul(&t, &a, &c);
    mpq_mul_2exp(t.val, t.val, 2);
    rat_sub(&d, &d, &t);
    q.complex = rat_is_negative(&d);
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

static void quad_roots_clear(QuadRoots* q) {
    rat_clear(&q->p); rat_clear(&q->r); mpz_clear(q->m);
}

// "p + " / "p - " / "" / "-" ahead of the radical term
static void append_root_prefix(StrBuf* sb, const QuadRoots* q, const char* sign) {
    if (mpq_sgn(q->p.val) != 0) {
        sb_append_owned(sb, rat_to_string(&q->p));
        sb_appendf(sb, " %s ", sign);
    } else if (strcmp(sign, "+") != 0) {
        sb_append(sb, strcmp(sign, "-") == 0 ? "-" : "+/- ");
    }
}

static void append_quadratic(StrBuf* sb, const IrrationalFactor* e, bool* first) {
    QuadRoots q = quad_roots(e->factor);
    if (q.complex) {
        separator(sb, first);
        sb_append(sb, "x = ");
        append_root_prefix(sb, &q, "+/-");
        if (!rat_is_one(&q.r) || mpz_cmp_ui(q.m, 1) != 0) { // "i" alone when the term is 1
            append_radical(sb, &q.r, q.m);
            sb_append(sb, "*");
        }
        sb_append(sb, "i");
        append_multiplicity(sb, e->multiplicity);
    } else {
        const char* signs[2] = { "-", "+" }; // intervals are ascending: smaller root first
        for (int i = 0; i < e->num_real; i++) {
            separator(sb, first);
            sb_append(sb, "x = ");
            append_root_prefix(sb, &q, signs[i]);
            append_radical(sb, &q.r, q.m);
            append_approx(sb, &e->real[i]);
            append_multiplicity(sb, e->multiplicity);
        }
    }
    quad_roots_clear(&q);
}

static void append_general(StrBuf* sb, const IrrationalFactor* e, bool* first) {
    for (int i = 0; i < e->num_real; i++) {
        separator(sb, first);
        sb_append(sb, "x");
        append_approx(sb, &e->real[i]);
        append_multiplicity(sb, e->multiplicity);
    }
    int nonreal = polynomial_degree(e->factor) - e->num_real;
    if (nonreal == 0) return;
    separator(sb, first);
    sb_appendf(sb, "%d non-real roots of (", nonreal);
    sb_append_owned(sb, poly_to_string(e->factor));
    sb_append(sb, ")");
    append_multiplicity(sb, e->multiplicity);
}

char* solutions_to_string(const Solutions* sol) {
    StrBuf sb = sb_new();
    bool first = true;
    for (int i = 0; i < sol->num_rational; i++) {
        separator(&sb, &first);
        sb_append(&sb, "x = ");
        sb_append_owned(&sb, rat_to_string(&sol->rational_roots[i]));
    }
    for (int i = 0; i < sol->num_irr; i++) {
        if (polynomial_degree(sol->irr[i].factor) == 2) append_quadratic(&sb, &sol->irr[i], &first);
        else append_general(&sb, &sol->irr[i], &first);
    }
    if (first) sb_append(&sb, "No solutions found");
    return sb_finish(&sb);
}

void free_solutions(Solutions* sol) {
    if (!sol) return;
    for (int i = 0; i < sol->num_rational; i++) rat_clear(&sol->rational_roots[i]);
    free(sol->rational_roots);
    for (int i = 0; i < sol->num_irr; i++) {
        RootIntervalList L = { sol->irr[i].real, sol->irr[i].num_real };
        root_interval_list_free(&L);
        free_polynomial(sol->irr[i].factor);
    }
    free(sol->irr);
    free(sol->irrational_roots);
    for (int i = 0; i < sol->num_irreducible; i++) free_polynomial(sol->irreducible_factors[i]);
    free(sol->irreducible_factors);
    free(sol);
}
