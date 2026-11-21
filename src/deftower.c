// Poles and tan crossings of f = N/D in θ = exp(v) or θ = tan(H) on a finite range (deftower.h).
// Points a + b·π are compared exactly when they agree in b, else with π to 50 digits; points
// with atan(r) or ln(r) (r ≠ 0, 1) are located in double precision with a margin, and a point
// inside the margin of an end is reported as undecided.
#include "deftower.h"
#include "kpoly.h"
#include "roots.h"
#include <math.h>
#include <stdlib.h>

#define MAX_POINTS 100000 // crossings of tan(H) on one range
#define MARGIN 1e-9       // relative distance from an end below which a pole is undecided

void pi_50(Rational* p) { rat_parse_literal(p, "3.14159265358979323846264338327950288419716939937511"); }

// --- points a + b·π ---

typedef struct {
    Rational a, b;
} PiPt;

static void pp_init(PiPt* p) { rat_init(&p->a); rat_init(&p->b); }
static void pp_clear(PiPt* p) { rat_clear(&p->a); rat_clear(&p->b); }

static void pp_of_end(PiPt* p, const Endpoint* e) {
    pp_init(p);
    rat_set(e->pi ? &p->b : &p->a, e->x);
}

static double pp_double(const PiPt* p) { return rat_to_double(&p->a) + rat_to_double(&p->b) * M_PI; }

// sign(p - q); 2 when π to 50 digits cannot decide it (|p - q| < 1e-50·|Δb|)
static int pp_cmp(const PiPt* p, const PiPt* q) {
    Rational da, db, pi, eps, v1, v2;
    rat_init(&da); rat_init(&db); rat_init(&pi); rat_init(&eps); rat_init(&v1); rat_init(&v2);
    rat_sub(&da, &p->a, &q->a);
    rat_sub(&db, &p->b, &q->b);
    int s;
    if (rat_is_zero(&db)) s = rat_is_zero(&da) ? 0 : rat_is_negative(&da) ? -1 : 1;
    else {                                                // da + db·[π - 1e-50, π + 1e-50]
        pi_50(&pi);
        rat_parse_literal(&eps, "1e-50");
        rat_sub(&v1, &pi, &eps); rat_mul(&v1, &v1, &db); rat_add(&v1, &v1, &da);
        rat_add(&v2, &pi, &eps); rat_mul(&v2, &v2, &db); rat_add(&v2, &v2, &da);
        bool pos = !rat_is_negative(&v1) && !rat_is_zero(&v1) && !rat_is_negative(&v2) && !rat_is_zero(&v2);
        bool neg = rat_is_negative(&v1) && rat_is_negative(&v2);
        s = pos ? 1 : neg ? -1 : 2;
    }
    rat_clear(&da); rat_clear(&db); rat_clear(&pi); rat_clear(&eps); rat_clear(&v1); rat_clear(&v2);
    return s;
}

// where p lies: -1 outside, 0 at an end, 1 strictly inside, 2 undecided
static int place(const PiPt* p, const PiPt ends[2]) {
    int l = pp_cmp(p, &ends[0]), h = pp_cmp(p, &ends[1]);
    if (l == 2 || h == 2) return 2;
    if (l < 0 || h > 0) return -1;
    return l == 0 || h == 0 ? 0 : 1;
}

// --- the tower's data ---

// D's coefficients as a polynomial over Q (they must be constants), or NULL
static Polynomial* constant_poly(const LPoly* D) {
    Polynomial* p = create_polynomial();
    Rational c, d;
    rat_init(&c); rat_init(&d);
    bool ok = true;
    for (int i = 0; ok && i < D->n; i++) {
        const QFn* a = &D->a[i];
        ok = polynomial_degree(a->num) <= 0 && polynomial_degree(a->den) == 0;
        if (!ok || poly_is_zero(a->num)) continue;
        poly_coeff(&c, a->num, 0);
        poly_coeff(&d, a->den, 0);
        rat_div(&c, &c, &d);
        add_term(p, &c, D->lo + i);
    }
    rat_clear(&c); rat_clear(&d);
    if (!ok) { free_polynomial(p); return NULL; }
    return p;
}

// h0, h1 of a linear polynomial h1·x + h0; false otherwise
static bool linear(const Polynomial* v, Rational* h0, Rational* h1) {
    if (polynomial_degree(v) != 1) return false;
    poly_coeff(h0, v, 0);
    poly_coeff(h1, v, 1);
    return true;
}

// Σ_j a_j(x0)·r^j for r = 0 or 1 (N at a rational point)
static bool numerator_zero_at(const LPoly* N, const Rational* x0, int r) {
    Rational acc, v, d;
    rat_init(&acc); rat_init(&v); rat_init(&d);
    for (int i = 0; i < N->n; i++) {
        if (r == 0 && N->lo + i != 0) continue;
        poly_eval_rat(&v, N->a[i].num, x0);
        poly_eval_rat(&d, N->a[i].den, x0);
        rat_div(&v, &v, &d);
        rat_add(&acc, &acc, &v);
    }
    bool z = rat_is_zero(&acc);
    rat_clear(&acc); rat_clear(&v); rat_clear(&d);
    return z;
}

// P / x^m with x ∤ the result
static Polynomial* without_zero_root(const Polynomial* P) {
    Polynomial* q = poly_copy(P);
    Rational z;
    rat_init(&z);
    for (poly_coeff(&z, q, 0); rat_is_zero(&z) && polynomial_degree(q) > 0; poly_coeff(&z, q, 0)) {
        Polynomial* x = create_polynomial();
        Rational one = rat_one_val();
        add_term(x, &one, 1);
        Polynomial* d = poly_divide(q, x);
        free_polynomial(q); free_polynomial(x);
        q = d;
        rat_clear(&one);
    }
    rat_clear(&z);
    return q;
}

// the real roots of P (square-free), each in an interval of width <= 1e-15
static RootIntervalList fine_roots(const Polynomial* P) {
    Polynomial *dp = poly_derivative(P), *g = poly_gcd(P, dp), *sf = poly_divide(P, g);
    RootIntervalList L = isolate_real_roots(sf);
    Rational w;
    rat_init(&w);
    rat_parse_literal(&w, "1e-15");
    for (int i = 0; i < L.count; i++) refine_root(&L.items[i], sf, &w);
    rat_clear(&w);
    free_polynomial(dp); free_polynomial(g); free_polynomial(sf);
    return L;
}

#define DIVERGES "the integrand has a pole in [a, b], so the integral diverges"
#define UNDECIDED "a pole of the integrand lies too close to a bound to be decided"
#define REMOVABLE "the integrand's denominator vanishes in [a, b] where its numerator may too; this is not supported yet"

// a numeric point interval [x1, x2] against the ends: the verdict for a pole there
static const char* numeric_pole(double x1, double x2, const double end[2]) {
    double m = MARGIN * (1 + fabs(x1) + fabs(x2));
    if (x2 < end[0] - m || x1 > end[1] + m) return NULL;
    return x1 > end[0] + m && x2 < end[1] - m ? DIVERGES : UNDECIDED;
}

// --- tan(H) ---

typedef struct {
    Rational h0, h1;
    PiPt ends[2];
    double end[2];
} TanCtx;

typedef struct {
    long k1, k2;
} KRange;

// the k range covering [end0, end1] for H = h1·x + h0 at H = (k + off)·π
static KRange k_range(const TanCtx* T, double off) {
    double h1 = rat_to_double(&T->h1), h0 = rat_to_double(&T->h0);
    double a = (h1 * T->end[0] + h0) / M_PI - off, b = (h1 * T->end[1] + h0) / M_PI - off;
    return (KRange){ (long)floor(fmin(a, b)) - 1, (long)ceil(fmax(a, b)) + 1 };
}

// point k of the family off, x = α + (k + off)/h1·π with α = -h0/h1 (exact a + b·π)
static PiPt tan_point(const TanCtx* T, long k, double off) {
    PiPt p;
    pp_init(&p);
    rat_div(&p.a, &T->h0, &T->h1);
    rat_neg(&p.a, &p.a);
    Rational kb = rat_from_si(2 * k + (off > 0 ? 1 : 0), 2);
    rat_div(&p.b, &kb, &T->h1);
    rat_clear(&kb);
    return p;
}

// θ = ∞ points: poles when deg N > deg D, else the crossings (recorded)
static const char* tan_infinities(const TanCtx* T, const ExpForm* F, TowerPoints* out) {
    KRange K = k_range(T, 0.5);
    if (K.k2 - K.k1 > MAX_POINTS) return "too many crossings of tan in [a, b]";
    bool pole = kp_deg(&F->f) > kp_deg(&F->den);
    const char* why = NULL;
    for (long k = K.k1; k <= K.k2 && !why; k++) {
        PiPt p = tan_point(T, k, 0.5);
        int w = place(&p, T->ends);
        if (w == 2) why = UNDECIDED;
        else if (w >= 0 && pole) why = DIVERGES;
        else if (w == 1) { out->cross = realloc(out->cross, (size_t)(out->n + 1) * sizeof(double)); out->cross[out->n++] = pp_double(&p); }
        else if (w == 0) { if (pp_cmp(&p, &T->ends[0]) == 0) out->at_lo = true; else out->at_hi = true; }
        pp_clear(&p);
    }
    return why;
}

// θ = 0 (D(0) = 0): x = α + k/h1·π; the rational point α is a pole unless N vanishes there
static const char* tan_zero_root(const TanCtx* T, const ExpForm* F) {
    KRange K = k_range(T, 0);
    if (K.k2 - K.k1 > MAX_POINTS) return "too many crossings of tan in [a, b]";
    const char* why = NULL;
    for (long k = K.k1; k <= K.k2 && !why; k++) {
        PiPt p = tan_point(T, k, 0);
        int w = place(&p, T->ends);
        if (w == 2) why = UNDECIDED;
        else if (w >= 0) why = k == 0 && numerator_zero_at(&F->f, &p.a, 0) ? REMOVABLE : DIVERGES;
        pp_clear(&p);
    }
    return why;
}

// a root r ≠ 0 of D in [r1, r2]: x = (atan(r) + kπ - h0)/h1, transcendental, a pole
static const char* tan_root(const TanCtx* T, double r1, double r2) {
    double h1 = rat_to_double(&T->h1), h0 = rat_to_double(&T->h0);
    KRange K = k_range(T, 0);
    const char* why = NULL;
    for (long k = K.k1 - 1; k <= K.k2 + 1 && !why; k++) {
        double a = (atan(r1) + k * M_PI - h0) / h1, b = (atan(r2) + k * M_PI - h0) / h1;
        why = numeric_pole(fmin(a, b), fmax(a, b), T->end);
    }
    return why;
}

// f's tower form, D's constant coefficients as a polynomial over Q, and the ends
typedef struct {
    const ExpForm* F;
    const Polynomial* D;
    const TowerEnds* E;
} TowerIn;

static const char* tan_points(const TowerIn* in, TowerPoints* out) {
    const ExpForm* F = in->F;
    const Polynomial* D = in->D;
    const TowerEnds* E = in->E;
    TanCtx T;
    rat_init(&T.h0); rat_init(&T.h1);
    const char* why = linear(F->t.v, &T.h0, &T.h1) ? NULL : "tan of a nonlinear argument in a denominator is not supported yet";
    pp_of_end(&T.ends[0], E->lo); pp_of_end(&T.ends[1], E->hi);
    T.end[0] = pp_double(&T.ends[0]); T.end[1] = pp_double(&T.ends[1]);
    if (!why) why = tan_infinities(&T, F, out);
    Rational z;
    rat_init(&z);
    poly_coeff(&z, D, 0);
    if (!why && rat_is_zero(&z)) why = tan_zero_root(&T, F);
    Polynomial* Dn = without_zero_root(D);                // its roots are nonzero
    RootIntervalList L = why ? (RootIntervalList){ NULL, 0 } : fine_roots(Dn);
    for (int i = 0; i < L.count && !why; i++) why = tan_root(&T, rat_to_double(&L.items[i].lo), rat_to_double(&L.items[i].hi));
    root_interval_list_free(&L);
    free_polynomial(Dn);
    rat_clear(&z); rat_clear(&T.h0); rat_clear(&T.h1);
    pp_clear(&T.ends[0]); pp_clear(&T.ends[1]);
    return why;
}

// --- exp(v) ---

// positive roots of D: r = 1 at the rational point -v0/v1, others at (ln r - v0)/v1
static const char* exp_points(const TowerIn* in) {
    const ExpForm* F = in->F;
    const Polynomial* D = in->D;
    const TowerEnds* E = in->E;
    RootIntervalList L = fine_roots(D);
    Rational v0, v1, one = rat_one_val(), d;
    rat_init(&v0); rat_init(&v1); rat_init(&d);
    bool lin = linear(F->t.v, &v0, &v1);
    PiPt ends[2];
    pp_of_end(&ends[0], E->lo); pp_of_end(&ends[1], E->hi);
    double end[2] = { pp_double(&ends[0]), pp_double(&ends[1]) };
    poly_eval_rat(&d, D, &one);
    bool one_root = rat_is_zero(&d);                      // θ = 1 is a root
    const char* why = NULL;
    for (int i = 0; i < L.count && !why; i++) {
        const RootInterval* r = &L.items[i];
        bool positive = !rat_is_negative(&r->hi) && !rat_is_zero(&r->hi);
        if (!positive) continue;                          // exp(v) > 0 never meets it
        if (!lin) why = "exp of a nonlinear argument in a denominator with positive roots is not supported yet";
        else if (one_root && rat_cmp(&r->lo, &one) <= 0 && rat_cmp(&r->hi, &one) >= 0) {
            PiPt p;                                       // θ = 1: x = -v0/v1
            pp_init(&p);
            rat_div(&p.a, &v0, &v1);
            rat_neg(&p.a, &p.a);
            int w = place(&p, ends);
            if (w >= 0) why = numerator_zero_at(&F->f, &p.a, 1) ? REMOVABLE : DIVERGES;
            pp_clear(&p);
        } else {
            double a = (log(fmax(rat_to_double(&r->lo), 1e-300)) - rat_to_double(&v0)) / rat_to_double(&v1);
            double b = (log(rat_to_double(&r->hi)) - rat_to_double(&v0)) / rat_to_double(&v1);
            why = numeric_pole(fmin(a, b), fmax(a, b), end);
        }
    }
    root_interval_list_free(&L);
    rat_clear(&v0); rat_clear(&v1); rat_clear(&one); rat_clear(&d);
    pp_clear(&ends[0]); pp_clear(&ends[1]);
    return why;
}

const char* tower_points(const ExpForm* F, const TowerEnds* ends, TowerPoints* out) {
    *out = (TowerPoints){ NULL, 0, false, false };
    Polynomial* D = constant_poly(&F->den);
    if (!D) return "a denominator in exp or tan with coefficients depending on x is not supported yet";
    TowerIn in = { F, D, ends };
    const char* why = F->t.kind == TOWER_TAN ? tan_points(&in, out) : exp_points(&in);
    free_polynomial(D);
    return why;
}

void tower_points_free(TowerPoints* p) { free(p->cross); p->cross = NULL; p->n = 0; }
