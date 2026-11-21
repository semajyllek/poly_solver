// The logarithmic part of ∫ r/(d·y), y = sqrt(U), deg U = 2, for an irreducible monic d coprime to U.
// Each root ρ of d gives c·∫ dx/((x - ρ)·y), c = r(ρ)/d'(ρ), and for any w with w² = U(ρ)
//   ∫ dx/((x - ρ)·y) = -(1/w)·ln((2U(ρ) + U'(ρ)·(x - ρ) + 2w·y)/(x - ρ)).
// With t = -c/(2w), the term is 2t·ln G; s = t² = c²/(4U(ρ)) lies in Q(ρ), and
//   s rational: w = c/(2√s) for every ρ, so the terms combine through the norm
//               E = Π (L(ρ) + B(ρ)·ν), ν = y/√s, into one ln (s > 0) or one atan (s < 0);
//   s generating Q(ρ): Σ over all roots t of R(t²) (R the minimal polynomial of s) of
//               t·ln(L(ρ) + 2w·y), the ln(x - ρ) parts cancelling between t and -t.
#include "sqrtlog.h"
#include "factor.h"
#include "matrix.h"
#include "residue.h"
#include <stdlib.h>

// --- arithmetic in Q(ρ) = Q[x]/d ---

// a·b mod d
static Polynomial* mulmod(const Polynomial* a, const Polynomial* b, const Polynomial* d) {
    Polynomial *p = poly_multiply(a, b), *r = poly_mod(p, d);
    free_polynomial(p);
    return r;
}

// 1/a mod d (a ≢ 0, d irreducible)
static Polynomial* invmod(const Polynomial* a, const Polynomial* d) {
    PolyBezout bz = poly_xgcd(a, d);
    Polynomial* s = poly_mod(bz.s, d);
    free_bezout(&bz);
    return s;
}

// p mod d in column j of M (the coordinates in 1, ρ, ..., ρ^(n-1))
static void put_column(Matrix* M, int j, const Polynomial* p) {
    for (int i = 0; i < M->rows; i++) poly_coeff(&MAT_AT(M, i, j), p, i);
}

// Π over the roots of d of e(ρ): the determinant of multiplication by e
static void norm_of(Rational* out, const Polynomial* e, const Polynomial* d) {
    int n = polynomial_degree(d);
    Matrix* M = matrix_create(n, n);
    Polynomial* p = poly_mod(e, d);
    for (int j = 0; j < n; j++) {
        put_column(M, j, p);
        Polynomial *x = create_polynomial(), *q;
        Rational one = rat_one_val();
        add_term(x, &one, 1);
        rat_clear(&one);
        q = mulmod(p, x, d);
        free_polynomial(x); free_polynomial(p);
        p = q;
    }
    matrix_det(out, M);
    free_polynomial(p);
    matrix_free(M);
}

// the polynomial Σ v[i]·x^i
static Polynomial* from_vector(const Rational* v, int n) {
    Polynomial* p = create_polynomial();
    for (int i = n - 1; i >= 0; i--)
        if (!rat_is_zero(&v[i])) add_term(p, &v[i], i);
    return p;
}

// --- s generating Q(ρ): elements as polynomials in s ---

// the basis change to powers of s: Pinv·coords(e) gives e = Σ a_i·s^i (NULL when s does not generate)
typedef struct {
    const Polynomial *d, *s;
    Matrix* Pinv;
    int n;
} PowerBasis;

static bool power_basis(PowerBasis* B) {
    Matrix* P = matrix_create(B->n, B->n);
    Rational one = rat_one_val();
    Polynomial* p = poly_constant(&one);
    rat_clear(&one);
    for (int j = 0; j < B->n; j++) {
        put_column(P, j, p);
        Polynomial* q = mulmod(p, B->s, B->d);
        free_polynomial(p);
        p = q;
    }
    free_polynomial(p);
    B->Pinv = matrix_inverse(P);
    matrix_free(P);
    return B->Pinv != NULL;
}

// e ∈ Q(ρ) as a polynomial in s of degree < n
static Polynomial* in_s(const PowerBasis* B, const Polynomial* e) {
    Polynomial* em = poly_mod(e, B->d);
    Rational* v = malloc((size_t)B->n * sizeof(Rational));
    Rational c, t;
    rat_init(&c); rat_init(&t);
    for (int i = 0; i < B->n; i++) {
        rat_init(&v[i]);
        for (int k = 0; k < B->n; k++) {
            poly_coeff(&c, em, k);
            rat_mul(&t, &MAT_AT(B->Pinv, i, k), &c);
            rat_add(&v[i], &v[i], &t);
        }
    }
    Polynomial* p = from_vector(v, B->n);
    for (int i = 0; i < B->n; i++) rat_clear(&v[i]);
    free(v);
    free_polynomial(em);
    rat_clear(&c); rat_clear(&t);
    return p;
}

// p(t²) mod q
static Polynomial* at_square(const Polynomial* p, const Polynomial* q) {
    Polynomial* r = create_polynomial();
    for (const Term* t = p->head; t; t = t->next) add_term(r, &t->coeff, 2 * t->exponent);
    Polynomial* m = poly_mod(r, q);
    free_polynomial(r);
    return m;
}

// the minimal polynomial R of s (monic, degree n): s^n = Σ a_i·s^i
static Polynomial* min_poly(const PowerBasis* B) {
    Polynomial *sn = poly_pow(B->s, (unsigned)B->n), *a = in_s(B, sn);
    Polynomial* R = create_polynomial();
    Rational one = rat_one_val();
    add_term(R, &one, B->n);
    rat_clear(&one);
    Polynomial* r = poly_subtract(R, a);
    free_polynomial(sn); free_polynomial(a); free_polynomial(R);
    return r;
}

// the pieces of L(ρ) + 2w·y = c0 + c1·x + c2·y as polynomials in ρ (the caller reduces mod d)
typedef struct {
    Polynomial *c0, *c1, *B;   // c0 = 2U(ρ) - ρ·U'(ρ), c1 = U'(ρ), B = c = r(ρ)/d'(ρ)
} LogPieces;

static LogPieces log_pieces(const SqrtLog* in) {
    Polynomial *Ur = poly_mod(in->U, in->d), *dU = poly_derivative(in->U), *dUr = poly_mod(dU, in->d);
    Polynomial *x = create_polynomial(), *xdU, *twoU, *dd = poly_derivative(in->d), *idd = invmod(dd, in->d);
    Rational one = rat_one_val(), two = rat_from_int(2);
    add_term(x, &one, 1);
    xdU = mulmod(x, dUr, in->d);
    twoU = poly_scale(Ur, &two);
    LogPieces P = { poly_subtract(twoU, xdU), dUr, mulmod(in->r, idd, in->d) };
    free_polynomial(Ur); free_polynomial(dU); free_polynomial(x); free_polynomial(xdU);
    free_polynomial(twoU); free_polynomial(dd); free_polynomial(idd);
    rat_clear(&one); rat_clear(&two);
    return P;
}

static void log_pieces_free(LogPieces* P) { free_polynomial(P->c0); free_polynomial(P->c1); free_polynomial(P->B); }

// k[0] + k[1]·x + k[2]·y as an expression (NULL when all are zero)
static Expr* linear_xy(const Rational* k, const Polynomial* U) {
    Polynomial* p = create_polynomial();
    if (!rat_is_zero(&k[1])) add_term(p, &k[1], 1);
    if (!rat_is_zero(&k[0])) add_term(p, &k[0], 0);
    Expr* e = poly_is_zero(p) ? NULL : create_poly_expr(p);
    if (!e) free_polynomial(p);
    if (rat_is_zero(&k[2])) return e;
    Expr* y = create_mul(create_const_rat(&k[2]), create_sqrt(create_poly_expr(poly_copy(U))));
    return e ? create_add(e, y) : y;
}

// one RootSum over an irreducible factor q of R(t²): S(t) = c0 + c1·x + 2w·y, w = -t·B/(2s), s = t²
static Expr* rootsum_over(const Polynomial* q, const LogPieces* Ps, const SqrtLog* in) {
    Polynomial *t2 = create_polynomial(), *t1 = create_polynomial();
    Rational one = rat_one_val(), mone, k[3];
    rat_init(&mone); rat_init(&k[0]); rat_init(&k[1]); rat_init(&k[2]);
    rat_neg(&mone, &one);
    add_term(t2, &one, 2); add_term(t1, &mone, 1);           // t², -t
    Polynomial *c0 = at_square(Ps->c0, q), *c1 = at_square(Ps->c1, q), *B = at_square(Ps->B, q);
    Polynomial *is = invmod(t2, q), *Bt = mulmod(B, t1, q), *w2 = mulmod(Bt, is, q);   // 2w = -t·B/s
    int n = polynomial_degree(q);
    Expr** c = calloc((size_t)n, sizeof(Expr*));
    for (int i = 0; i < n; i++) {
        poly_coeff(&k[0], c0, i); poly_coeff(&k[1], c1, i); poly_coeff(&k[2], w2, i);
        c[i] = linear_xy(k, in->U);
    }
    Polynomial* qm = poly_copy(q);
    poly_make_monic(qm);
    free_polynomial(t2); free_polynomial(t1); free_polynomial(c0); free_polynomial(c1); free_polynomial(B);
    free_polynomial(is); free_polynomial(Bt); free_polynomial(w2);
    rat_clear(&one); rat_clear(&mone); rat_clear(&k[0]); rat_clear(&k[1]); rat_clear(&k[2]);
    return create_rootsum_terms(qm, c);
}

// s generating Q(ρ): Σ over the irreducible factors q of R(t²) of RootSum(q, t·ln S(t))
static Expr* generic_part(const LogPieces* Ps, const PowerBasis* PB, const SqrtLog* in) {
    LogPieces Pt = { in_s(PB, Ps->c0), in_s(PB, Ps->c1), in_s(PB, Ps->B) };   // in s, then s = t²
    Polynomial *R = min_poly(PB), *R2 = create_polynomial();
    for (const Term* t = R->head; t; t = t->next) add_term(R2, &t->coeff, 2 * t->exponent);
    Factorization* f = factorize(R2);
    Expr* e = NULL;
    for (int i = 0; i < f->count; i++) e = expr_append_sum(e, rootsum_over(f->factors[i], &Pt, in));
    free_factorization(f);
    free_polynomial(R); free_polynomial(R2);
    log_pieces_free(&Pt);
    return e;
}

// --- s rational: one norm ---

// the polynomial through (i, ys[i]) for i = 0..n (Lagrange)
static Polynomial* interpolate(const Rational* ys, int n) {
    Polynomial* p = create_polynomial();
    for (int i = 0; i <= n; i++) {
        Rational one = rat_one_val(), den = rat_one_val(), t;
        rat_init(&t);
        Polynomial* l = poly_constant(&one);
        for (int k = 0; k <= n; k++) {
            if (k == i) continue;
            Polynomial* f = create_polynomial();
            Rational mk = rat_from_int(-k), ik = rat_from_int(i - k);
            add_term(f, &one, 1); add_term(f, &mk, 0);
            Polynomial* g = poly_multiply(l, f);
            free_polynomial(l); free_polynomial(f);
            l = g;
            rat_mul(&den, &den, &ik);
            rat_clear(&mk); rat_clear(&ik);
        }
        rat_div(&t, &ys[i], &den);
        Polynomial *lt = poly_scale(l, &t), *sum = poly_add(p, lt);
        free_polynomial(p); free_polynomial(l); free_polynomial(lt);
        p = sum;
        rat_clear(&one); rat_clear(&den); rat_clear(&t);
    }
    return p;
}

// E(x, ν) = Π_ρ (c0 + c1·x + B·ν) as Σ_l e[l](x)·ν^l (l <= n), by interpolation on a grid
static Polynomial** norm_in_nu(const LogPieces* Ps, const Polynomial* d) {
    int n = polynomial_degree(d);
    Rational* v = malloc((size_t)(n + 1) * sizeof(Rational));
    Polynomial** Ej = malloc((size_t)(n + 1) * sizeof(Polynomial*));       // E(x, j) in x
    for (int i = 0; i <= n; i++) rat_init(&v[i]);
    for (int j = 0; j <= n; j++) {
        for (int i = 0; i <= n; i++) {
            Rational xi = rat_from_int(i), nj = rat_from_int(j);
            Polynomial *a = poly_scale(Ps->c1, &xi), *b = poly_scale(Ps->B, &nj), *s1 = poly_add(Ps->c0, a), *s2 = poly_add(s1, b);
            norm_of(&v[i], s2, d);
            free_polynomial(a); free_polynomial(b); free_polynomial(s1); free_polynomial(s2);
            rat_clear(&xi); rat_clear(&nj);
        }
        Ej[j] = interpolate(v, n);
    }
    Polynomial** e = malloc((size_t)(n + 1) * sizeof(Polynomial*));
    for (int l = 0; l <= n; l++) e[l] = create_polynomial();
    for (int j = 0; j <= n; j++) {                            // Σ_j E(x, j)·ℓ_j(ν)
        for (int i = 0; i <= n; i++) mpq_set_si(v[i].val, i == j, 1);   // the unit vector e_j
        Polynomial* lj = interpolate(v, n);
        for (int l = 0; l <= n; l++) {
            Rational c;
            rat_init(&c);
            poly_coeff(&c, lj, l);
            Polynomial *t = poly_scale(Ej[j], &c), *s = poly_add(e[l], t);
            free_polynomial(t); free_polynomial(e[l]);
            e[l] = s;
            rat_clear(&c);
        }
        free_polynomial(lj); free_polynomial(Ej[j]);
    }
    for (int i = 0; i <= n; i++) rat_clear(&v[i]);
    free(v); free(Ej);
    return e;
}

// Σ e[l]·ν^l, l <= n
typedef struct {
    Polynomial** e;
    int n;
} NuPoly;

// E = P0 + P1·ν with ν² = U/s
static void split_norm(const NuPoly* E, const Polynomial* Us, Polynomial** P) {
    Polynomial** e = E->e;
    int n = E->n;
    P[0] = create_polynomial(); P[1] = create_polynomial();
    Rational one = rat_one_val();
    Polynomial* pw = poly_constant(&one);                    // (U/s)^(l/2)
    rat_clear(&one);
    for (int l = 0; l <= n; l++) {
        Polynomial *t = poly_multiply(e[l], pw), *s = poly_add(P[l % 2], t);
        free_polynomial(t); free_polynomial(P[l % 2]);
        P[l % 2] = s;
        if (l % 2) { Polynomial* q = poly_multiply(pw, Us); free_polynomial(pw); pw = q; }
    }
    free_polynomial(pw);
}

// p·y, or y alone for p = 1
static Expr* times_y(const Polynomial* p, const Polynomial* U) {
    Expr* y = create_sqrt(create_poly_expr(poly_copy(U)));
    Rational one = rat_one_val(), c;
    rat_init(&c);
    poly_coeff(&c, p, 0);
    bool unit = polynomial_degree(p) == 0 && rat_eq(&c, &one);
    rat_clear(&one); rat_clear(&c);
    return unit ? y : create_mul(create_poly_expr(poly_copy(p)), y);
}

// P = g·Pp with Pp primitive and a positive leading coefficient
static Polynomial* primitive(const Polynomial* P, Rational* g) {
    Rational c = poly_content(P), l;
    rat_init(&l);
    poly_coeff(&l, P, polynomial_degree(P));
    rat_abs(g, &c);
    if (rat_is_negative(&l)) rat_neg(g, g);
    rat_clear(&c); rat_clear(&l);
    Rational ig;
    rat_init(&ig);
    mpq_inv(ig.val, g->val);
    Polynomial* p = poly_scale(P, &ig);
    rat_clear(&ig);
    return p;
}

// s > 0, √s = k·√m: -2√s·(ln|P0 + P1·y/√s| - ln|d|), the ln argument divided by the content of P0;
// s < 0, √|s| = k·√m: 2√|s|·atan(√|s|·P0/(P1·y)) with P0, P1 primitive
static Expr* norm_terms(Polynomial** P, const Rational* s, const SqrtLog* in) {
    Rational as, k, c, g0, g1, t;
    rat_init(&as); rat_init(&k); rat_init(&c); rat_init(&g0); rat_init(&g1); rat_init(&t);
    rat_abs(&as, s);
    mpz_t m;
    mpz_init(m);
    rat_sqrt_split(&as, &k, m);
    Polynomial *p0 = poly_is_zero(P[0]) ? NULL : primitive(P[0], &g0), *p1 = primitive(P[1], &g1);
    if (!p0) rat_set(&g0, &g1);
    Expr* e;
    if (!rat_is_negative(s)) {
        mpq_set_z(c.val, m);
        rat_mul(&c, &c, &k);
        rat_mul(&c, &c, &g0);
        rat_div(&c, &g1, &c);                                 // (g1/g0)/√s = c·√m
        Expr* yt = surd_times(&c, m, times_y(p1, in->U));
        Expr* arg = p0 ? create_add(create_poly_expr(poly_copy(p0)), yt) : yt;
        rat_add(&t, &k, &k);
        rat_neg(&t, &t);
        e = surd_times(&t, m, create_ln(arg));
        rat_neg(&t, &t);
        e = expr_append_sum(e, surd_times(&t, m, create_ln(create_poly_expr(poly_copy(in->d)))));
    } else {
        rat_mul(&c, &k, &g0);
        rat_div(&c, &c, &g1);
        Rational p00;
        rat_init(&p00);
        if (p0) poly_coeff(&p00, p0, 0);
        bool unit = p0 && polynomial_degree(p0) == 0 && rat_is_one(&p00);           // c·√m/(P1·y)
        rat_clear(&p00);
        Expr* q = unit ? create_div(surd_times(&c, m, NULL), times_y(p1, in->U))
                       : surd_times(&c, m, create_div(p0 ? create_poly_expr(poly_copy(p0)) : create_const(0), times_y(p1, in->U)));
        rat_add(&t, &k, &k);
        e = surd_times(&t, m, create_atan(q));
    }
    if (p0) free_polynomial(p0);
    free_polynomial(p1);
    mpz_clear(m);
    rat_clear(&as); rat_clear(&k); rat_clear(&c); rat_clear(&g0); rat_clear(&g1); rat_clear(&t);
    return e;
}

static Expr* rational_part(const LogPieces* Ps, const Rational* s, const SqrtLog* in) {
    int n = polynomial_degree(in->d);
    Polynomial** e = norm_in_nu(Ps, in->d);
    Rational is;
    rat_init(&is);
    mpq_inv(is.val, s->val);
    Polynomial *Us = poly_scale(in->U, &is), *P[2];
    split_norm(&(NuPoly){ e, n }, Us, P);
    Expr* r = poly_is_zero(P[1]) ? create_const(0) : norm_terms(P, s, in);   // P1 = 0: no y, the part is 0
    for (int l = 0; l <= n; l++) free_polynomial(e[l]);
    free(e);
    free_polynomial(Us); free_polynomial(P[0]); free_polynomial(P[1]);
    rat_clear(&is);
    return r;
}

Expr* sqrt_log_part(const SqrtLog* in, const char** why) {
    LogPieces Ps = log_pieces(in);
    Polynomial *B2 = mulmod(Ps.B, Ps.B, in->d), *U4 = poly_mod(in->U, in->d);
    Rational four = rat_from_int(4);
    Polynomial *U4s = poly_scale(U4, &four), *iU = invmod(U4s, in->d), *s = mulmod(B2, iU, in->d);   // s = c²/(4U(ρ))
    Expr* e = NULL;
    PowerBasis PB = { in->d, s, NULL, polynomial_degree(in->d) };
    if (polynomial_degree(s) <= 0) {
        Rational s0;
        rat_init(&s0);
        poly_coeff(&s0, s, 0);
        e = rational_part(&Ps, &s0, in);
        rat_clear(&s0);
    } else if (power_basis(&PB)) e = generic_part(&Ps, &PB, in);
    else *why = "the residues of this square-root integrand generate only part of the field of the denominator's roots (not supported)";
    if (PB.Pinv) matrix_free(PB.Pinv);
    log_pieces_free(&Ps);
    free_polynomial(B2); free_polynomial(U4); free_polynomial(U4s); free_polynomial(iU); free_polynomial(s);
    rat_clear(&four);
    return e;
}
