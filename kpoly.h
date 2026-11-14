#ifndef KPOLY_H
#define KPOLY_H

#include "expfield.h"

// Polynomials in θ over K = Q(x), stored as LPoly with lo >= 0: the Euclidean ring K[θ] that
// Hermite reduction and the residue test of the Risch algorithm work in. These are the
// algebraic operations (θ as an indeterminate); the derivation is lpoly_derivative.

typedef struct {
    LPoly q, r;
} KDiv;

// s·a + t·b = c with deg s < deg b
typedef struct {
    LPoly s, t;
} KBezout;

// d = Π f[i]^(i+1), f[i] monic, square-free and pairwise coprime (a part may be 1)
typedef struct {
    LPoly* f;
    int m;
} KSquareFree;

// n/d in lowest terms with d monic in θ (and, for θ = exp(v), θ ∤ d: θ powers live in n)
typedef struct {
    LPoly n, d;
} KFrac;

int kp_deg(const LPoly* a);                            // -1 for 0
LPoly kp_one(void);
bool kp_is_one(const LPoly* a);
LPoly kp_copy(const LPoly* a);
LPoly kp_shift(const LPoly* a, int k);                 // a·θ^k
LPoly kp_scale(const LPoly* a, const QFn* c);          // c·a
LPoly kp_monic(const LPoly* a);                        // a / lc(a), a ≠ 0
LPoly kp_dtheta(const LPoly* a);                       // d/dθ
KDiv kp_divmod(const LPoly* a, const LPoly* b);       // b ≠ 0
LPoly kp_exact_div(const LPoly* a, const LPoly* b);   // a/b when b | a
LPoly kp_gcd(const LPoly* a, const LPoly* b);         // monic (0 only when both are 0)
KBezout kp_bezout(const LPoly* a, const LPoly* b, const LPoly* c); // gcd(a, b) = 1, b ≠ 0
KSquareFree kp_square_free(const LPoly* d);           // Yun, d monic and nonconstant
void kdiv_free(KDiv* d);
void kbezout_free(KBezout* z);
void ksf_free(KSquareFree* s);

KFrac kfrac_make(LPoly n, LPoly d, TowerKind kind);   // takes both; d ≠ 0
KFrac kfrac_add(const KFrac* a, const KFrac* b, TowerKind kind);
KFrac kfrac_sub(const KFrac* a, const KFrac* b, TowerKind kind);
KFrac kfrac_derivative(const KFrac* f, const Tower* t);
void kfrac_free(KFrac* f);

LPoly kp_pow(const LPoly* a, int k);                   // a^k, k >= 0
LPoly kp_primitive(const LPoly* a);                    // a·c, c ∈ Q(x), coefficients coprime in Z[x], lc > 0
Expr* kfrac_expr(const KFrac* f, const Tower* t);     // as (num)/(den) over Q[x], den from the top power

#endif
