#include "solver.h"
#include "factor.h"
#include "rootform.h"
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

#define TINY_ROOT 1e-4 // below this, roots are refined relatively and printed in scientific form

static bool is_tiny(const Rational* q) { return fabs(rat_to_double(q)) < TINY_ROOT; }

// a root near 0 (x^3 - 10^30·x + 1 has one near 10^-30): bisect on until the interval keeps
// one sign and its width is 10^-12 of its size, so the digits shown are significant
static void refine_relative(RootInterval* r, const Polynomial* f) {
    Rational w, lo, hi;
    rat_init(&w); rat_init(&lo); rat_init(&hi);
    for (int step = 0; step < 4000; step++) {                 // a nonzero root: always ends far sooner
        rat_abs(&lo, &r->lo); rat_abs(&hi, &r->hi);
        rat_sub(&w, &r->hi, &r->lo);
        bool one_sign = mpq_sgn(r->lo.val) * mpq_sgn(r->hi.val) > 0;
        Rational* near = rat_cmp(&lo, &hi) < 0 ? &lo : &hi;
        mpq_div_2exp(near->val, near->val, 40);                // 2^-40 ≈ 10^-12 of the root
        if (one_sign && rat_cmp(&w, near) <= 0) break;
        mpq_div_2exp(w.val, w.val, 1);                         // one more bisection
        refine_root(r, f, &w);
    }
    rat_clear(&w); rat_clear(&lo); rat_clear(&hi);
}

// isolate and refine the real roots of an irreducible factor of degree >= 2
static void solve_irrational(Solutions* sol, const Polynomial* f, int mult) {
    RootIntervalList L = isolate_real_roots(f);
    Rational width; // 10^-(digits+2): well below the last digit shown
    rat_init(&width);
    mpz_ui_pow_ui(mpq_denref(width.val), 10, DECIMAL_DIGITS + 2);
    mpz_set_ui(mpq_numref(width.val), 1);
    for (int i = 0; i < L.count; i++) {
        refine_root(&L.items[i], f, &width);
        if (is_tiny(&L.items[i].lo) || is_tiny(&L.items[i].hi)) refine_relative(&L.items[i], f);
    }
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
    if (is_tiny(&mid)) {                                      // 10 significant digits: 1.000000000e-30
        mpf_t m;
        mpf_init2(m, 256);
        mpf_set_q(m, mid.val);
        char* s = NULL;
        gmp_asprintf(&s, "%.9Fe", m);
        sb_append(sb, s);
        free(s);
        mpf_clear(m);
    } else sb_append_owned(sb, rat_to_decimal(&mid, DECIMAL_DIGITS));
    rat_clear(&mid);
}

static void append_multiplicity(StrBuf* sb, int mult) {
    if (mult > 1) sb_appendf(sb, " (multiplicity %d)", mult);
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
    QuadRoots q = poly_quad_roots(e->factor);
    if (q.is_complex) {
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

// an exact root part as text, parenthesized when it is a sum
static void append_part(StrBuf* sb, const Expr* e) {
    bool sum = e->type == EXPR_BINOP && (e->binop.op == OP_ADD || e->binop.op == OP_SUB);
    if (sum) sb_append(sb, "(");
    sb_append_owned(sb, expr_to_string(e));
    if (sum) sb_append(sb, ")");
}

// position of θ (units of π) among the roots: 0 for -ρ (θ = 1), 1 for +ρ (θ = 0),
// 2 for the upper member of a conjugate pair (0 < θ < 1), 3 for the lower member
static int root_class(const Rational* theta) {
    Rational one = rat_one_val();
    int c = rat_cmp(theta, &one);
    rat_clear(&one);
    if (c == 0) return 0;
    if (rat_is_zero(theta)) return 1;
    return c < 0 ? 2 : 3;
}

// "x = 2^(1/3) ~ 1.2599210498" with the approximation from the certified interval
static void append_known_real(StrBuf* sb, const RootSet* rs, const Rational* theta) {
    Expr* x = root_part(rs, theta, TRIG_COS);
    sb_append(sb, "x = ");
    append_part(sb, x);
    free_expr(x);
}

// "x = a +/- b * i ~ (1.0911 +/- 0.6300i)" for the pair ρ·e^(±iπθ); the decimals are display only
static void append_known_pair(StrBuf* sb, const RootSet* rs, const Rational* theta) {
    Expr *re = root_part(rs, theta, TRIG_COS), *im = root_part(rs, theta, TRIG_SIN);
    bool ok = true;
    double a = re ? expr_eval(re, 0, &ok) : 0, b = expr_eval(im, 0, &ok);
    sb_append(sb, "x = ");
    if (re) { append_part(sb, re); sb_append(sb, " +/- "); }
    else sb_append(sb, "+/- ");
    if (!(im->type == EXPR_CONST && rat_is_one(&im->constant))) { append_part(sb, im); sb_append(sb, " * "); }
    sb_append(sb, "i");
    if (fabs(a) < TINY_ROOT && fabs(b) < TINY_ROOT) sb_appendf(sb, " ~ (%.9e +/- %.9ei)", a, b); // significant digits
    else sb_appendf(sb, " ~ (%.10f +/- %.10fi)", a, b);
    if (re) free_expr(re);
    free_expr(im);
}

// every root of a binomial or rescaled-cyclotomic factor exactly; false for other factors
static bool append_known(StrBuf* sb, const IrrationalFactor* e, bool* first) {
    RootSet rs;
    if (!known_roots(e->factor, &rs)) return false;
    int real = 0;
    for (int cls = 0; cls <= 2; cls++)                  // -ρ, +ρ (ascending), then the pairs
        for (int i = 0; i < rs.count; i++) {
            Rational theta = rat_from_si(rs.js[i], rs.L);
            if (root_class(&theta) == cls) {
                separator(sb, first);
                if (cls < 2) { append_known_real(sb, &rs, &theta); append_approx(sb, &e->real[real++]); }
                else append_known_pair(sb, &rs, &theta);
                append_multiplicity(sb, e->multiplicity);
            }
            rat_clear(&theta);
        }
    root_set_free(&rs);
    return true;
}

static void append_general(StrBuf* sb, const IrrationalFactor* e, bool* first) {
    if (append_known(sb, e, first)) return;
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
