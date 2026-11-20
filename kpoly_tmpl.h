// Polynomial algorithms over a coefficient field, written once and instantiated per type at
// compile time (no void* dispatch on the hot paths): include this file after defining
//   P            the polynomial type (lo >= 0 powers, coefficient array a[n])
//   C            the coefficient type
//   CTX          a context argument type passed to every operation (the tower kind for KFrac)
//   PFX(name)    the prefix for the generated function names
//   P_DEG(p)               degree, -1 for 0
//   P_LC(p)                const C*: the leading coefficient
//   P_TERM(x, c, k)        c·θ^k for const C* c (copied)
//   P_ADD(x, a, b, sign)   a ± b
//   P_MUL(x, a, b)         a·b
//   P_ONE(x), P_ZERO(x)    1 and 0
//   P_FREE(p)
//   P_COEF(p, k)           const C*: the coefficient of θ^k, k within [lo, lo + n)
//   P_LO(p), P_N(p)
//   C_MUL(x, a, b), C_INV(x, a), C_FREE(a), C_SCALE_INT(x, a, k) (k·a)
// and KDIV {P q, r}, KBEZ {P s, t}, KSF {P* f; int m} as the result types. The generated
// functions mirror kpoly.h: copy, scale, monic, dtheta, divmod, exact_div, gcd (monic),
// inverse, bezout and square_free.

static inline P PFX(copy)(CTX x, const P* a) {
    (void)x;                                                   // unused by some instances
    P z = P_ZERO(x), r = P_ADD(x, a, &z, 1);
    P_FREE(&z);
    return r;
}

// c·a
static inline P PFX(scale)(CTX x, const P* a, const C* c) {
    (void)x;
    P t = P_TERM(x, c, 0), r = P_MUL(x, a, &t);
    P_FREE(&t);
    return r;
}

// a / lc(a), a ≠ 0
static inline P PFX(monic)(CTX x, const P* a) {
    (void)x;
    C inv = C_INV(x, P_LC(a));
    P r = PFX(scale)(x, a, &inv);
    C_FREE(&inv);
    return r;
}

// d/dθ
static inline P PFX(dtheta)(CTX x, const P* a) {
    (void)x;
    P r = P_ZERO(x);
    for (int i = 0; i < P_N(a); i++) {
        int k = P_LO(a) + i;
        if (k == 0) continue;
        C c = C_SCALE_INT(x, P_COEF(a, k), k);
        P t = P_TERM(x, &c, k - 1), s = P_ADD(x, &r, &t, 1);
        C_FREE(&c);
        P_FREE(&r); P_FREE(&t);
        r = s;
    }
    return r;
}

// a = q·b + r with deg r < deg b, from the top one term at a time
static inline KDIV PFX(divmod)(CTX x, const P* a, const P* b) {
    KDIV d = { P_ZERO(x), PFX(copy)(x, a) };
    C ib = C_INV(x, P_LC(b));
    while (P_DEG(&d.r) >= P_DEG(b)) {
        C c = C_MUL(x, P_LC(&d.r), &ib);
        P t = P_TERM(x, &c, P_DEG(&d.r) - P_DEG(b));          // the next quotient term
        C_FREE(&c);
        P tb = P_MUL(x, &t, b), r = P_ADD(x, &d.r, &tb, -1), q = P_ADD(x, &d.q, &t, 1);
        P_FREE(&d.r); P_FREE(&d.q); P_FREE(&t); P_FREE(&tb);
        d.r = r; d.q = q;
    }
    C_FREE(&ib);
    return d;
}

static inline P PFX(exact_div)(CTX x, const P* a, const P* b) {
    KDIV d = PFX(divmod)(x, a, b);
    P_FREE(&d.r);
    return d.q;
}

// monic gcd (0 only when both are 0)
static inline P PFX(gcd)(CTX x, const P* a, const P* b) {
    P r0 = PFX(copy)(x, a), r1 = PFX(copy)(x, b);
    while (P_N(&r1)) {                                         // (r0, r1) <- (r1, r0 mod r1)
        KDIV d = PFX(divmod)(x, &r0, &r1);
        P_FREE(&d.q); P_FREE(&r0);
        r0 = r1; r1 = d.r;
    }
    P g = P_N(&r0) ? PFX(monic)(x, &r0) : r0;
    if (P_N(&r0)) P_FREE(&r0);
    P_FREE(&r1);                                               // the zero remainder still owns its array
    return g;
}

// s with s·a ≡ 1 (mod b), by the extended Euclidean algorithm (only the s column is needed)
static inline P PFX(inverse)(CTX x, const P* a, const P* b) {
    P r0 = PFX(copy)(x, a), r1 = PFX(copy)(x, b), s0 = P_ONE(x), s1 = P_ZERO(x);
    while (P_N(&r1)) {
        KDIV d = PFX(divmod)(x, &r0, &r1);
        P qs = P_MUL(x, &d.q, &s1), s2 = P_ADD(x, &s0, &qs, -1);
        P_FREE(&qs); P_FREE(&d.q); P_FREE(&r0); P_FREE(&s0);
        r0 = r1; r1 = d.r;
        s0 = s1; s1 = s2;
    }
    C inv = C_INV(x, P_LC(&r0));                               // r0 is the constant gcd
    P s = PFX(scale)(x, &s0, &inv);
    C_FREE(&inv);
    P_FREE(&r0); P_FREE(&r1); P_FREE(&s0); P_FREE(&s1);
    return s;
}

// s·a + t·b = c with deg s < deg b, for gcd(a, b) = 1 and b ≠ 0
static inline KBEZ PFX(bezout)(CTX x, const P* a, const P* b, const P* c) {
    P inv = PFX(inverse)(x, a, b), cs = P_MUL(x, c, &inv);
    KDIV sd = PFX(divmod)(x, &cs, b);                          // s = c·a⁻¹ mod b
    P sa = P_MUL(x, &sd.r, a), rest = P_ADD(x, c, &sa, -1);
    KBEZ z = { sd.r, PFX(exact_div)(x, &rest, b) };            // t = (c - s·a)/b
    P_FREE(&inv); P_FREE(&cs); P_FREE(&sd.q); P_FREE(&sa); P_FREE(&rest);
    return z;
}

// Yun: with c = gcd(d, d_θ), w = d/c, y = d_θ/c, each step peels off the part of multiplicity i
static inline KSF PFX(square_free)(CTX x, const P* d) {
    KSF s = { malloc((size_t)(P_DEG(d) + 1) * sizeof(P)), 0 };
    P dp = PFX(dtheta)(x, d), c = PFX(gcd)(x, d, &dp);
    P w = PFX(exact_div)(x, d, &c), y = PFX(exact_div)(x, &dp, &c);
    P_FREE(&dp); P_FREE(&c);
    while (P_DEG(&w) > 0) {
        P wp = PFX(dtheta)(x, &w), z = P_ADD(x, &y, &wp, -1);
        P g = PFX(gcd)(x, &w, &z);
        P nw = PFX(exact_div)(x, &w, &g), ny = PFX(exact_div)(x, &z, &g);
        s.f[s.m++] = g;                                        // the part of multiplicity m
        P_FREE(&wp); P_FREE(&z); P_FREE(&w); P_FREE(&y);
        w = nw; y = ny;
    }
    P_FREE(&w); P_FREE(&y);
    return s;
}
