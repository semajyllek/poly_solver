# polysolver Internals

A comprehensive guide to how every algorithm, data structure, and module in polysolver works, illustrated with toy examples.

---

## Table of Contents

1. [System Overview](#1-system-overview)
2. [Module Dependency Graph](#2-module-dependency-graph)
3. [Memory Ownership Model](#3-memory-ownership-model)
4. [Rational Arithmetic](#4-rational-arithmetic)
5. [Polynomial Algebra](#5-polynomial-algebra)
6. [Expression Trees](#6-expression-trees)
7. [The Shunting-Yard Parser](#7-the-shunting-yard-parser)
8. [Canonicalization](#8-canonicalization)
9. [Polynomial Factoring](#9-polynomial-factoring)
10. [Equation Solving](#10-equation-solving)
11. [Symbolic Calculus](#11-symbolic-calculus)
12. [CLI Pipeline](#12-cli-pipeline)

---

## 1. System Overview

polysolver is a symbolic algebra engine written in C. It takes a string like `"x^2 - 5x + 6"`, parses it into an internal tree representation, reduces it to a canonical polynomial form, and then applies operations—simplification, factoring, differentiation, integration, or equation solving.

The data flows through a pipeline:

```
Input string
    │
    ▼
┌──────────┐
│  Lexer   │  token.c — breaks string into tokens (numbers, variables, operators)
└──────────┘
    │  stream of Token*
    ▼
┌──────────────────┐
│ Shunting-Yard    │  parser.c — converts token stream to expression tree
│ Parser           │
└──────────────────┘
    │  Expr* (tree of BINOP, VAR, CONST, NEG nodes)
    ▼
┌──────────────────┐
│ Canonicalization  │  expression.c — reduces tree to EXPR_POLY or EXPR_RATIONAL_FN
└──────────────────┘
    │  Expr* (containing Polynomial* or num/den Polynomial* pair)
    ▼
┌──────────────────┐
│ Operation        │  factor.c / solver.c / expression.c
│ (factor, solve,  │  — operates on the canonical polynomial form
│  diff, integrate)│
└──────────────────┘
    │
    ▼
  Output string
```

---

## 2. Module Dependency Graph

```
strbuf.c/h                ← growable strings; every printer uses it (no fixed buffers)
rational_num.c/h          ← GMP mpq_t wrapper, exact literal parsing, decimal output
    │
    ▼
polynomial.c/h            ← sparse Q[x]: arithmetic, gcd, xgcd, pow, exact evaluation
    │
    ├──► poly_zp.c/h      ← dense F_p[x]: arithmetic, distinct/equal-degree factoring
    ├──► poly_zx.c/h      ← dense Z[x] (mpz_t): exact division, primitive part, mod p^k
    │        │
    │        ▼
    │    zx_factor.c      ← Berlekamp–Zassenhaus: prime choice, Hensel lifting, recombination
    │    modgcd.c/h       ← gcd and inverse over Q by CRT + rational reconstruction, checked exactly
    │        │
    ▼        ▼
factor.c/h                ← content, Yun square-free decomposition, full factorization
    │
    ├──► roots.c/h        ← Descartes/VCA real-root isolation, exact bisection refinement
    │        │
    │        ▼
    ├──► solver.c/h       ← exact rational roots, radicals for quadratics, certified intervals
    │
expression.c/h            ← Expr tree, canonicalization, printing, numeric evaluation
    │
    ├──► integrate.c      ← Hermite reduction + per-factor log/atan terms
    │        └──► logpart.c, subres.c ← Lazard–Rioboo–Trager for factors of degree >= 3
    ├──► partial_frac.c/h ← partial fractions over Q-irreducible factors (uses matrix)
    │
token.c, parser.c         ← lexer (exact literal text) and shunting-yard parser
matrix.c/h                ← exact Gaussian elimination (rational_num, strbuf)
cli/main.c                ← CLI wiring (the library sources are in src/)
```

---

## 3. Memory Ownership Model

The entire codebase follows one rule:

> **Caller owns inputs. Every function returns new allocations. No function mutates or frees its inputs.**

This means:

```c
Polynomial* a = create_polynomial();
add_term(a, rat_from_int(3), 2);  // 3x^2

Polynomial* b = create_polynomial();
add_term(b, rat_from_int(1), 2);  // x^2

Polynomial* sum = poly_add(a, b); // returns NEW polynomial 4x^2
// a and b are untouched — caller must free all three
free_polynomial(a);
free_polynomial(b);
free_polynomial(sum);
```

Ownership transfer happens in only one place: the `extract_poly` and `extract_ratfn` helpers inside canonicalization, which null out pointers in the source `Expr` to transfer ownership of the inner `Polynomial*` without copying.

The `Rational` type wraps a GMP `mpq_t`. It must be initialized (`rat_init`, or a constructor like `rat_from_int`) and cleared (`rat_clear`); arithmetic writes into a pre-initialized result. `rat_to_string()` and `rat_to_decimal()` return heap strings the caller frees.

---

## 4. Rational Arithmetic

### 4.1 Why exact rationals?

Floating-point arithmetic accumulates rounding errors. Consider computing `1/3 + 1/3 + 1/3`:

```
With doubles:  0.333333... + 0.333333... + 0.333333... = 0.999999...  (≠ 1)
With rationals: 1/3 + 1/3 + 1/3 = 3/3 = 1  (exact)
```

Exact arithmetic is essential for polynomial GCD (where a single rounding error can change the degree of the result) and for the rational root theorem (where we need to know if a polynomial evaluates to *exactly* zero at a candidate root).

### 4.2 The `Rational` struct

```c
typedef struct {
    mpq_t val;   // GMP arbitrary-precision rational, always canonical
} Rational;
```

GMP keeps every value in lowest terms with a positive denominator, so equality is structural and there is no overflow: `rat_add`, `rat_mul` and friends are thin wrappers over `mpq_add`, `mpq_mul`, etc.

### 4.3 Exact literals in, exact text out

Numbers never pass through `double`. The lexer finds a literal's extent with `rat_parse_literal`, which accepts `digits[.digits][(e|E)[+-]digits]` and builds the value exactly:

```
"1.25"   → 125 · 10^-2 = 5/4
"2e-2"   → 2 · 10^-2   = 1/50
"10000000000000000000000" → stays exact (it used to saturate to 2^63 - 1)
"2e"     → 2, then the variable e   (an exponent needs digits)
```

Matrix entries use the same parser, so `[[1/3, 0.5]]` is exact too. On output, `rat_to_decimal(r, d)` prints a truncated decimal with `d` places using only integer arithmetic, and all printers append into a growable `StrBuf`, so coefficients of any size print in full.

### 4.4 Square roots of rationals

`rat_sqrt_split(q, &r, m)` writes `sqrt(q) = r · sqrt(m)` with `m` an integer stripped of square factors (trial division up to 10^4, then a perfect-square test on what is left). The solver and integrator use it to print radicals such as `1/2*sqrt(3)`.

### 4.5 Exponentiation by squaring

`rat_pow_int` uses the binary exponentiation algorithm (also called "exponentiation by squaring"). To compute `base^exp`:

```
(2/3)^5:
    exp = 5 = 101 in binary
    
    result = 1/1, base = 2/3
    
    bit 0 (1): result = 1 * 2/3 = 2/3,     base = 2/3 * 2/3 = 4/9
    bit 1 (0): result unchanged,             base = 4/9 * 4/9 = 16/81
    bit 2 (1): result = 2/3 * 16/81 = 32/243
    
    Answer: 32/243
```

This computes the result in O(log n) multiplications instead of O(n).

### 4.6 Comparison

Two normalized rationals `a/b` and `c/d` (where `b, d > 0`) are compared by cross-multiplying: compare `a*d` versus `c*b`. Since denominators are positive, this preserves the ordering.

---

## 5. Polynomial Algebra

### 5.1 Data structure: sorted sparse linked list

A polynomial is represented as a singly-linked list of `Term` nodes, sorted in **descending order by exponent**:

```c
typedef struct Term {
    Rational coeff;    // coefficient (exact rational)
    int exponent;      // degree of this term
    int var_id;        // 0 = "x" (placeholder for future multivariate)
    struct Term* next; // pointer to next term (lower degree)
} Term;

typedef struct Polynomial {
    Term* head;        // highest-degree term (or NULL for zero polynomial)
} Polynomial;
```

This is a **sparse** representation: only terms with nonzero coefficients are stored. The polynomial `x^100 + 1` has just two nodes, not 101.

**Example:** The polynomial `3x^2 - 4x + 5` is stored as:

```
Polynomial { head ──► }

    ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
    │ coeff: 3/1  │     │ coeff: -4/1 │     │ coeff: 5/1  │
    │ exp:   2    │────►│ exp:   1    │────►│ exp:   0    │────► NULL
    └─────────────┘     └─────────────┘     └─────────────┘
```

The zero polynomial is simply `Polynomial { head: NULL }`.

### 5.2 The `add_term` insertion algorithm

`add_term` is the fundamental building block. It inserts a term into a polynomial while maintaining sorted order and combining like terms. Here's a step-by-step trace:

**Example:** Adding `2x` to the polynomial `3x^2 + 5`:

```
Starting state:
    head ──► [3, exp=2] ──► [5, exp=0] ──► NULL

add_term(poly, 2/1, 1):
    1. Create new_term: [2, exp=1]
    2. Walk the list to find position:
       - current = [3, exp=2], exp=2 > 1 → advance
         previous = [3, exp=2], current = [5, exp=0]
       - current = [5, exp=0], exp=0 < 1 → stop (insert here)
    3. No existing term with exp=1, so insert between:
       previous→next = new_term
       new_term→next = current

Result:
    head ──► [3, exp=2] ──► [2, exp=1] ──► [5, exp=0] ──► NULL
```

**Example:** Adding `-3x^2` to the above (like-term combination):

```
add_term(poly, -3/1, 2):
    1. Walk the list:
       - current = [3, exp=2], exp=2 == 2 → match found!
    2. Add coefficients: 3 + (-3) = 0
    3. Coefficient is zero → remove the node entirely

Result:
    head ──► [2, exp=1] ──► [5, exp=0] ──► NULL
```

This automatic cancellation means the list never contains zero-coefficient terms.

### 5.3 Polynomial addition

`add_term` scans from the head, so building a polynomial with it costs O(n²). The core operations avoid it: they build their result with an appender (a tail pointer) in descending order, O(1) per term. `poly_add` and `poly_subtract` merge the two sorted lists, combining equal exponents and dropping zeros:

```
(3x^2 + 4x) + (5x^2 + 1):
    exp 2 in both → 3 + 5 = 8   → [8, exp=2]
    exp 1 only in the first      → [8, exp=2] → [4, exp=1]
    exp 0 only in the second     → [8, exp=2] → [4, exp=1] → [1, exp=0]

Result: 8x^2 + 4x + 1
```

`poly_copy`, `poly_scale`, `poly_derivative` and `poly_integral` append in the same way.

### 5.4 Polynomial multiplication

Every term of poly1 times every term of poly2, accumulated in a dense window of `mpq_t` over the exponents `lo1 + lo2 .. deg1 + deg2`, then read back from the top:

```
(2x + 3) * (x + 1):
    window[2] += 2·1   window[1] += 2·1   window[1] += 3·1   window[0] += 3·1

Result: 2x^2 + 5x + 3
```

Time complexity: O(n*m) coefficient operations, with no list walks. Very sparse factors (a window more than 64·n·m + 1024 wide, like `x^100000·y`) insert term by term instead.

### 5.5 Polynomial long division (`poly_divmod`)

This implements the classical polynomial long division algorithm, returning both quotient and remainder. The algorithm repeatedly divides the leading term of the remainder by the leading term of the divisor:

**Example:** `(x^3 + 2x^2 + 3x + 4) ÷ (x + 1)`:

```
dividend: x^3 + 2x^2 + 3x + 4
divisor:  x + 1

Step 1: leading terms: x^3 / x = x^2
    quotient so far: x^2
    subtract x^2 * (x + 1) = x^3 + x^2 from remainder:
    remainder: (x^3 + 2x^2 + 3x + 4) - (x^3 + x^2) = x^2 + 3x + 4

Step 2: leading terms: x^2 / x = x
    quotient so far: x^2 + x
    subtract x * (x + 1) = x^2 + x from remainder:
    remainder: (x^2 + 3x + 4) - (x^2 + x) = 2x + 4

Step 3: leading terms: 2x / x = 2
    quotient so far: x^2 + x + 2
    subtract 2 * (x + 1) = 2x + 2 from remainder:
    remainder: (2x + 4) - (2x + 2) = 2

Step 4: degree(remainder) = 0 < degree(divisor) = 1 → stop

Result: quotient = x^2 + x + 2, remainder = 2
```

The implementation works in place: the dividend goes into a dense window of `mpq_t`, and each step computes the quotient coefficient `top · lc(divisor)⁻¹` and subtracts it times the divisor's lower terms from the window, so no intermediate polynomial is built. The convenience wrappers `poly_divide` and `poly_mod` call `poly_divmod` and discard the unwanted part.

### 5.6 Polynomial GCD (Euclidean algorithm)

The polynomial GCD uses the same Euclidean algorithm as integer GCD, but with polynomial remainder instead of integer remainder:

```c
Polynomial* poly_gcd(poly1, poly2) {
    r0 = copy of poly1
    r1 = copy of poly2
    while r1 ≠ 0:
        r2 = r0 mod r1
        r0 = r1
        r1 = r2
    make r0 monic (divide by leading coefficient)
    return r0
}
```

The implementation also makes each remainder monic as it goes. That leaves the gcd unchanged and keeps the rationals smaller. `poly_gcd_fast` (§9.3) uses the modular algorithm instead once both inputs have degree 6 or more.

**Example:** `gcd(x^2 - 1, x^2 - 2x + 1)`:

We know `x^2 - 1 = (x-1)(x+1)` and `x^2 - 2x + 1 = (x-1)^2`, so the GCD should be `(x-1)`.

```
r0 = x^2 - 1
r1 = x^2 - 2x + 1

Iteration 1:
    r2 = (x^2 - 1) mod (x^2 - 2x + 1)
    Long division: (x^2 - 1) ÷ (x^2 - 2x + 1)
        leading: x^2/x^2 = 1
        subtract 1*(x^2 - 2x + 1) from (x^2 - 1):
        = (x^2 - 1) - (x^2 - 2x + 1) = 2x - 2
    r2 = 2x - 2
    r0 = x^2 - 2x + 1
    r1 = 2x - 2

Iteration 2:
    r2 = (x^2 - 2x + 1) mod (2x - 2)
    Long division: (x^2 - 2x + 1) ÷ (2x - 2)
        leading: x^2/(2x) = x/2
        subtract (x/2)*(2x-2) = x^2 - x:
        remainder = -x + 1
        leading: -x/(2x) = -1/2
        subtract (-1/2)*(2x-2) = -x + 1:
        remainder = 0
    r2 = 0
    r0 = 2x - 2
    r1 = 0

Loop ends. Make monic: (2x - 2) / 2 = x - 1.
Result: gcd = x - 1  ✓
```

Key design: the function **copies** its inputs before the Euclidean loop. The original `poly_gcd` had a bug where it freed `poly1` during the loop, destroying the caller's data.

### 5.7 Horner's method for evaluation

The naive approach to evaluating `3x^2 + 4x + 5` at `x = 2` computes each power independently: `3*(2^2) + 4*(2) + 5 = 12 + 8 + 5 = 25`. This requires computing `x^k` for each term.

Horner's method rewrites the polynomial as nested multiplications:

```
3x^2 + 4x + 5 = ((3)*x + 4)*x + 5
```

The algorithm walks the sorted linked list from highest to lowest degree. For each gap between consecutive exponents, it multiplies by `x`:

```
Evaluating 3x^2 + 4x + 5 at x = 2:

    Start: result = 3 (coefficient of x^2)
    
    Gap from exp=2 to exp=1: 1 step
        result = 3 * 2 = 6
    Add coefficient of x^1:
        result = 6 + 4 = 10
    
    Gap from exp=1 to exp=0: 1 step
        result = 10 * 2 = 20
    Add coefficient of x^0:
        result = 20 + 5 = 25
```

For sparse polynomials like `x^100 + 1`, Horner handles the gap efficiently—it multiplies by `x` 99 times between the two terms, then once more for the trailing gap to degree 0.

### 5.8 Calculus operations on polynomials

**Derivative** applies the power rule to each term: `d/dx(a*x^n) = a*n*x^(n-1)`. Constant terms (exponent = 0) are dropped.

```
d/dx(3x^2 + 4x + 5):
    3x^2 → 3*2 * x^1 = 6x
    4x   → 4*1 * x^0 = 4
    5    → dropped (exponent = 0)
Result: 6x + 4
```

**Integral** applies the reverse power rule: `∫a*x^n dx = (a/(n+1))*x^(n+1)`. The constant of integration is not tracked.

```
∫(6x + 4) dx:
    6x → (6/2)*x^2 = 3x^2
    4  → (4/1)*x^1 = 4x
Result: 3x^2 + 4x
```

The exact rational arithmetic is critical here: integrating `x^2` gives `(1/3)*x^3`, which is represented exactly as coefficient `{1, 3}`.

---

## 6. Expression Trees

### 6.1 The `Expr` tagged union

Before canonicalization, parsed expressions are represented as trees of `Expr` nodes. The `Expr` struct uses a **tagged union** — a type tag (`ExprType`) determines which member of the union is active:

```c
typedef struct Expr {
    ExprType type;       // discriminant tag
    union {
        Rational constant;          // EXPR_CONST: a number like 3 or 1/2
        char* var_name;             // EXPR_VAR: "x"
        struct {                    // EXPR_BINOP: binary operation
            BinOp op;               //   which operation (+, -, *, /, ^)
            struct Expr* left;      //   left operand
            struct Expr* right;     //   right operand
        } binop;
        struct Expr* operand;       // EXPR_NEG: unary negation
        Polynomial* poly;           // EXPR_POLY: canonical polynomial
        struct {                    // EXPR_RATIONAL_FN: canonical rational function
            Polynomial* num;
            Polynomial* den;
        } ratfn;
    };
} Expr;
```

### 6.2 Tree structure example

The expression `(x + 2) * x^3` is represented as:

```
        MUL
       /   \
     ADD    POW
    /   \   / \
  VAR  CONST VAR CONST
  "x"   2   "x"   3
```

In memory:

```
Expr{EXPR_BINOP, .binop = {
    op: OP_MUL,
    left:  Expr{EXPR_BINOP, .binop = {
        op: OP_ADD,
        left:  Expr{EXPR_VAR, .var_name = "x"},
        right: Expr{EXPR_CONST, .constant = {2, 1}}
    }},
    right: Expr{EXPR_BINOP, .binop = {
        op: OP_POW,
        left:  Expr{EXPR_VAR, .var_name = "x"},
        right: Expr{EXPR_CONST, .constant = {3, 1}}
    }}
}}
```

### 6.3 Unified binary operations

All five binary operations (ADD, SUB, MUL, DIV, POW) share the same struct layout—just a `BinOp` discriminant and two child pointers. This replaces an earlier design that had five separate structs (`add.left/add.right`, `mul.left/mul.right`, etc.), reducing boilerplate in every switch statement from 5 cases to 1.

### 6.4 Memory management via `free_expr`

`free_expr` recursively walks the tree by tag:

```c
free_expr(MUL node):
    → free_expr(left child)   // recurse
    → free_expr(right child)  // recurse
    → free(node itself)

free_expr(VAR node):
    → free(var_name string)   // strdup'd at creation
    → free(node itself)

free_expr(POLY node):
    → free_polynomial(poly)   // frees linked list
    → free(node itself)

free_expr(NULL):
    → no-op (safe to call on NULL)
```

---

## 7. The Shunting-Yard Parser

### 7.1 Why shunting-yard?

Dijkstra's shunting-yard algorithm converts an infix expression (like `2 + 3 * x`) into a structured form that respects operator precedence, in a single left-to-right pass. It uses two stacks: one for operands and one for operators.

The name comes from a railroad switchyard analogy: tokens arrive on one track, get routed through a switching area (the operator stack), and come out on another track in the correct order.

### 7.2 Data structures

The parser uses two type-safe stacks:

```c
// Operand stack — holds expression tree nodes
typedef struct { Expr** items; int top; int cap; } ExprStack;

// Operator stack — holds operator metadata
typedef struct {
    TokenType type;    // TOKEN_PLUS, TOKEN_MULTIPLY, etc.
    bool is_unary;     // distinguishes unary minus from binary minus
} OpEntry;
typedef struct { OpEntry* items; int top; int cap; } OpStack;
```

### 7.3 Precedence and associativity table

```
Level 3: power (^), unary minus   right-associative
Level 2: multiply (*), divide (/) left-associative
Level 1: add (+), subtract (-)    left-associative
Level 0: parentheses, sentinel
```

Power is right-associative so `x^2^3` parses as `x^(2^3)`, not `(x^2)^3`.

Unary minus shares the level of `^` and is also right-associative, so `-x^2` is `-(x^2)` and `x^-2` is `x^(-2)`, while `-2*3` still applies the sign first. An earlier version gave unary minus a higher level, so `-x^2` parsed as `(-x)^2`.

After parsing, `parse_and_simplify` rejects any unconsumed input (`x^2 = 4`, `2x $ 3`). A missing operand (`2 +`) or an unbalanced parenthesis is also an error, never a shortened expression.

### 7.4 The algorithm, step by step

Here's a detailed trace of parsing `2x^3 + 5x - 1`:

The lexer produces tokens: `2`, `x`, `^`, `3`, `+`, `5`, `x`, `-`, `1`

```
Processing token: 2 (NUMBER)
    prev_was_operand = false → no implicit mul
    Push CONST(2) onto expr stack
    prev_was_operand = true

    Expr stack:  [CONST(2)]
    Op stack:    []

Processing token: x (VARIABLE)
    prev_was_operand = true, x is operand-start → insert implicit *
        push_operator(*, prec=2): op stack empty → just push
    Push VAR("x") onto expr stack
    prev_was_operand = true

    Expr stack:  [CONST(2), VAR(x)]
    Op stack:    [* (prec 2)]

Processing token: ^ (POWER)
    prev_was_operand = true, ^ is not operand-start → no implicit mul
    push_operator(^, prec=3, right-assoc):
        top of op stack: * (prec 2), 2 < 3 → don't pop. Push ^.
    prev_was_operand = false

    Expr stack:  [CONST(2), VAR(x)]
    Op stack:    [*, ^]

Processing token: 3 (NUMBER)
    prev_was_operand = false → no implicit mul
    Push CONST(3) onto expr stack
    prev_was_operand = true

    Expr stack:  [CONST(2), VAR(x), CONST(3)]
    Op stack:    [*, ^]

Processing token: + (PLUS)
    prev_was_operand = true, + is not operand-start → no implicit mul
    push_operator(+, prec=1):
        top of op stack: ^ (prec 3), 3 > 1 → pop and apply:
            pop CONST(3), VAR(x) → create POW(VAR(x), CONST(3))
            push result
        top of op stack: * (prec 2), 2 > 1 → pop and apply:
            pop POW(x,3), CONST(2) → create MUL(CONST(2), POW(x,3))
            push result
        op stack empty → push +
    prev_was_operand = false

    Expr stack:  [MUL(2, x^3)]
    Op stack:    [+]

Processing token: 5 (NUMBER)
    prev_was_operand = false → no implicit mul
    Push CONST(5) onto expr stack
    prev_was_operand = true

    Expr stack:  [MUL(2, x^3), CONST(5)]
    Op stack:    [+]

Processing token: x (VARIABLE)
    prev_was_operand = true, x is operand-start → insert implicit *
        push_operator(*, prec=2):
            top: + (prec 1), 1 < 2 → don't pop. Push *.
    Push VAR("x") onto expr stack
    prev_was_operand = true

    Expr stack:  [MUL(2, x^3), CONST(5), VAR(x)]
    Op stack:    [+, *]

Processing token: - (MINUS)
    prev_was_operand = true → binary minus (not unary)
    push_operator(-, prec=1):
        top: * (prec 2), 2 > 1 → pop and apply:
            pop VAR(x), CONST(5) → create MUL(CONST(5), VAR(x))
            push result
        top: + (prec 1), 1 == 1 and left-assoc → pop and apply:
            pop MUL(5,x), MUL(2,x^3) → create ADD(MUL(2,x^3), MUL(5,x))
            push result
        op stack empty → push -

    Expr stack:  [ADD(2x^3, 5x)]
    Op stack:    [-]

Processing token: 1 (NUMBER)
    prev_was_operand = false → no implicit mul
    Push CONST(1) onto expr stack
    prev_was_operand = true

    Expr stack:  [ADD(2x^3, 5x), CONST(1)]
    Op stack:    [-]

End of input — pop remaining operators:
    pop - and apply:
        pop CONST(1), ADD(2x^3, 5x) → create SUB(ADD(2x^3, 5x), CONST(1))

Final expr stack: [SUB(ADD(2x^3, 5x), CONST(1))]
```

The result is the tree:

```
        SUB
       /    \
     ADD    CONST(1)
    /    \
  MUL    MUL
  / \    / \
 2  POW 5   x
    / \
   x   3
```

### 7.5 Implicit multiplication

The parser detects implicit multiplication by tracking whether the previous token was an operand. If it was, and the current token could start a new operand (number, variable, or left paren), a synthetic `TOKEN_MULTIPLY` is injected into the operator processing.

This handles cases like:

| Input | Interpretation |
|-------|---------------|
| `5x` | `5 * x` |
| `2x^3` | `2 * x^3` |
| `xy` | `x * y` |
| `(x+1)(x-1)` | `(x+1) * (x-1)` |
| `3(x+1)` | `3 * (x+1)` |

### 7.6 Unary minus detection

A minus token is **unary** if the previous token was not an operand — i.e., at the start of input, after an operator, or after an opening paren. Unary minus gets precedence level 4 (highest) and right-associativity.

```
Input: "-(x + 2)"

Token -: prev_was_operand = false → unary minus
    Push OpEntry{MINUS, is_unary=true} at precedence 4

Token (: push LEFT_PAREN on op stack

Token x: push VAR(x)

Token +: push_operator(+, prec=1)

Token 2: push CONST(2)

Token ): pop until LEFT_PAREN:
    pop + → apply: ADD(VAR(x), CONST(2))

End of input: pop remaining:
    pop unary - → apply: NEG(ADD(VAR(x), CONST(2)))
```

---

## 8. Canonicalization

### 8.1 What is canonicalization?

Canonicalization converts an arbitrary expression tree into one of two **canonical forms**:

1. **`EXPR_POLY`**: A single `Polynomial*` — the fully expanded, like-terms-combined form.
2. **`EXPR_RATIONAL_FN`**: A numerator and denominator, both `Polynomial*`, reduced by their GCD.

After canonicalization, `simplify()` is literally just `return canonicalize(expr)`.

### 8.2 The recursive algorithm

`canonicalize` processes the expression tree bottom-up:

**Base cases:**
- `EXPR_CONST` with value `c` → polynomial with single term: `c*x^0`
- `EXPR_VAR` → polynomial with single term: `1*x^1`

**Recursive cases:**
- `EXPR_NEG` → canonicalize the operand, then negate all coefficients
- `EXPR_BINOP` → canonicalize both children, then apply the operation on the resulting polynomials

### 8.3 Worked example

Canonicalizing `(x + 1)^2`:

```
Input tree:
    POW
    / \
  ADD  CONST(2)
  / \
VAR  CONST(1)
"x"

Step 1: Canonicalize left child (ADD):
    Canonicalize VAR("x") → poly: [1*x^1]
    Canonicalize CONST(1) → poly: [1*x^0]
    Apply OP_ADD:
        poly_add([1*x^1], [1*x^0]) → [1*x^1, 1*x^0]
    Result: EXPR_POLY representing (x + 1)

Step 2: Canonicalize right child (CONST(2)):
    → poly: [2*x^0]

Step 3: Apply OP_POW with integer exponent 2:
    lp = x + 1 (the polynomial)
    exp = 2

    result = copy of lp = x + 1
    iteration 1 (i=1, i < 2):
        result = poly_multiply(result, lp)
        = (x + 1) * (x + 1)
        = x^2 + 2x + 1

Result: EXPR_POLY representing x^2 + 2x + 1
```

### 8.4 Rational function handling

When a division appears, canonicalization produces an `EXPR_RATIONAL_FN`. The `wrap_ratfn` helper tries to simplify:

1. If the denominator is the constant polynomial `1`, return just the numerator as `EXPR_POLY`.
2. Otherwise, compute `gcd(numerator, denominator)` and divide both by it.
3. If the GCD-reduced denominator is now `1`, return as `EXPR_POLY`.

**Example:** Canonicalizing `(x^2 - 1) / (x - 1)`:

```
Numerator poly: x^2 - 1
Denominator poly: x - 1

wrap_ratfn:
    gcd(x^2 - 1, x - 1) = x - 1 (degree > 0, so reduce)
    numerator / gcd = (x^2 - 1) / (x - 1) = x + 1
    denominator / gcd = (x - 1) / (x - 1) = 1
    Denominator is now 1 → return EXPR_POLY: x + 1
```

### 8.5 Mixed rational function arithmetic

When one operand is `EXPR_POLY` and the other is `EXPR_RATIONAL_FN` (or both are rational functions), the `EXPR_POLY` is lifted to a rational function with denominator `1`, then standard fraction arithmetic is applied:

```
a/b + c/d = (a*d + c*b) / (b*d)
a/b * c/d = (a*c) / (b*d)
a/b ÷ c/d = (a*d) / (b*c)
```

The result is then passed through `wrap_ratfn` for GCD reduction.

---

## 9. Polynomial Factoring

### 9.1 Pipeline

```
f ∈ Q[x]
  │  Yun square-free decomposition (gcds over Q)
  ▼
f = lc · a_1 · a_2^2 · … · a_k^k        a_i square-free, pairwise coprime
  │  per a_i: primitive integer part, then Berlekamp–Zassenhaus
  ▼
irreducible primitive integer factors with multiplicities, sorted
  │  content = lc(f) / Π lc(factor)^mult
  ▼
f = content · Π factor^mult   (exactly)
```

The convention is the standard one: factors are primitive integer polynomials with positive leading coefficient, and all scaling lives in `content`. So `2x^2 - 3x + 1` prints as `(2x - 1)(x - 1)`, and `x^2/2 - 1/8` as `1/8(2x + 1)(2x - 1)`.

### 9.2 Yun's square-free decomposition

With `a0 = gcd(f, f')`, `b = f/a0`, `d = f'/a0 - b'`, each pass computes `a_i = gcd(b, d)`, then `b ← b/a_i` and `d ← (d/a_i) - b'`. Part `i` collects exactly the factors of multiplicity `i`. The gcds use `poly_gcd_fast` (§9.3):

```
f = (x - 1)(x^2 - 2)^2
parts: [x - 1, x^2 - 2]   → multiplicities 1 and 2
```

### 9.3 Modular gcd and inverse (`modgcd.c`)

Euclid over Q has coefficient blow-up: the intermediate remainders' numerators and denominators grow far beyond the answer's. For a dense degree-30 denominator, that made the gcds in Yun's algorithm and the inverse in `reduce_monic` take over 80% of the time. `modgcd.c` computes images in F_p for word-size primes, from `2³¹ − 1` downward, and combines them by CRT.
- **Unlucky primes.** For a gcd, only images of the least degree are kept. A prime that divides the resultant of the cofactors gives a longer image: it is dropped, or it restarts the accumulation when it came first. For an inverse, a prime where `a` is not a unit modulo `M` is skipped, as is any prime that divides a leading coefficient.
- **Reconstruction.** Rational reconstruction (half-extended Euclid, `|n|, d ≤ √(m/2)`) is tried each time the prime count has grown by a quarter, which keeps its total cost proportional to the last try.
- **Check.** When two tries agree, the candidate is checked exactly over Q. A gcd candidate must divide both inputs: a monic common divisor has degree ≤ deg gcd, and no image has smaller degree than the true gcd, so this proves equality. An inverse candidate must satisfy `s·a mod M = 1`.

A NULL result after 4000 primes falls back to Euclid, so results are always exact. `poly_gcd_fast` uses the modular route when both degrees are ≥ 6. Tests exercise both kinds of unlucky prime: `x` with `x − p₁p₂` (stable but wrong until the check), and a run of many primes in which the fifth is unlucky.

Measured CPU time, main → this, on dense random denominators (outputs byte-identical): degree 15: 0.12 → 0.03 s; degree 20: 0.65 → 0.12 s; degree 25: 3.9 → 0.45 s; degree 30: 14.3/10.1 → 1.1/0.66 s (13–15×).

### 9.4 Berlekamp–Zassenhaus (`zx_factor.c`)

For a primitive, square-free `f ∈ Z[x]`:

1. **Choose a prime.** Try odd primes until five are usable (p does not divide lc(f), and f stays square-free mod p). Factor f mod each one and keep the prime giving the fewest factors. One factor mod p proves f irreducible.
2. **Factor mod p** (`poly_zp.c`). Distinct-degree factorization isolates the product of all degree-d factors as `gcd(x^(p^d) - x, f)`. Cantor–Zassenhaus then splits each product with `gcd(a^((p^d-1)/2) - 1, g)` for pseudo-random `a`. The generator has a fixed seed, so output is deterministic.
3. **Hensel lift** to `P = p^k > 2·|lc|·(n+1)·2^n·max|a_i|`, a loosened Mignotte bound. The target is `f / lc(f) mod P`. Factors are split into two halves recursively; each pair `g·h` is lifted linearly with Bézout coefficients `s·g + t·h ≡ 1 (mod p)`:
   ```
   e  = (target - g·h) / M  mod p
   g += M · (e·t mod g)      h += M · (e·s mod h)      M ← M·p
   ```
4. **Recombine.** For subset sizes 1, 2, …, form `lc(f) · Π u_i mod P` in the symmetric range and take its primitive part. Accept it **only if it divides f exactly in Z[x]**. A lifting bug can therefore only under-factor; it can never return a wrong factor.

Examples that exercise recombination: `x^4 + 1` and `x^4 - 10x^2 + 1` are irreducible over Q but split modulo every prime. `x^8 + 98x^4 + 1` splits into two quartics.

---

## 10. Equation Solving

### 10.1 Strategy

`solve_polynomial` factors first. Linear factors give exact rational roots. Every irreducible factor of degree ≥ 2 becomes an `IrrationalFactor`, holding its certified real-root intervals.

### 10.2 Real-root isolation (`roots.c`)

The method is Vincent–Collins–Akritas bisection driven by Descartes' rule of signs. Everything is exact integer arithmetic:

- A power-of-two Cauchy bound `B` maps the positive roots of `f(Bt)` and of `f(-Bt)` into `(0, 1)`.
- The number of roots in `(0, 1)` is bounded by the sign variations of `(t+1)^n · g(1/(t+1))`, computed with a reversal and a Taylor shift. A count of 0 drops the interval, 1 isolates it, and more bisects via `2^n·g(t/2)` and its shift by 1.
- An exact root at a bisection point is recorded as a degenerate interval.

`refine_root` then bisects an isolating interval with exact evaluation (`poly_eval_rat`) until it is narrower than 10^-12. The printed decimal is the truncated midpoint. A root below 10^-4 in size (`x³ − 10³⁰·x + 1` has one near 10^-30) is bisected further, until the interval keeps one sign and is 2^-40 of the root wide, and is printed with 10 significant digits in scientific form (`1.000000000e-30`, through GMP's `mpf`); tiny complex roots print their decimals the same way.

### 10.3 Output

- Rational roots: `x = 3/2`.
- Quadratic factors print exact radicals: `x = 1/2 - 1/2*sqrt(5) ~ -0.6180339887`, or `x = -1/2 +/- 1/2*sqrt(3)*i` for complex roots.
- Binomial and rescaled-cyclotomic factors (`known_roots`, §11.4) print every root exactly: real ones as `x = 2^(1/3) ~ 1.2599210498`, conjugate pairs as `x = ρcosθ +/- ρsinθ * i ~ (a +/- bi)`, e.g. `x = cos(2*pi/5) +/- sin(2*pi/5) * i` for `x^5 - 1`. The real approximations come from the certified intervals; the complex ones are evaluated in double precision, for display only.
- Other higher-degree factors print their real roots as `x ~ -0.6823278038` and the number of non-real roots: `2 non-real roots of (x^3 + x + 1)`.
- Isolation can return an interval whose endpoint is a neighbouring exact root (`(0, 4)` beside `[0, 0]` for `x(x-1)(x+1)`). `refine_root` then takes the sign just right of `lo` from `f'(lo)`, which is nonzero because `f` is square-free.

For compatibility, `Solutions` still exposes `irrational_roots` (doubles from the certified intervals) and `irreducible_factors` (factors that have non-real roots).

---

## 11. Symbolic Calculus

### 11.1 Differentiation

`differentiate()` first canonicalizes the input expression, then dispatches:

**Polynomial case:** Delegates to `poly_derivative`, which applies the power rule term-by-term.

**Rational function case:** Applies the quotient rule:

```
d/dx(f/g) = (f'g - fg') / g^2
```

**Example:** Differentiating `(x^2 + 1) / (x - 1)`:

```
f  = x^2 + 1     f' = 2x
g  = x - 1       g' = 1

Numerator: f'g - fg' = 2x(x-1) - (x^2+1)(1)
         = 2x^2 - 2x - x^2 - 1
         = x^2 - 2x - 1

Denominator: g^2 = (x-1)^2 = x^2 - 2x + 1

Result: (x^2 - 2x - 1) / (x^2 - 2x + 1)
```

The result is then passed through `wrap_ratfn` which attempts GCD reduction. In this case, the GCD of the numerator and denominator is 1, so no simplification occurs.

### 11.2 Integration (`integrate.c`)

For `N/D`:

1. **Polynomial part.** Divide: `N = q·D + r`, and integrate `q` term by term.
2. **Rational part: Hermite reduction** (Mack's linear version; Bronstein, *Symbolic Integration I*, §2.2). This needs only gcds and the extended Euclidean algorithm, no factoring:
   ```
   Dm ← gcd(D, D'), Ds ← D/Dm
   while deg Dm > 0:
       Dm2 ← gcd(Dm, Dm'), Dms ← Dm/Dm2
       solve B·(-Ds·Dm'/Dm) + C·Dms = A
       A ← C - B'·Ds/Dms,  g ← g + B/Dm,  Dm ← Dm2
   ```
   The result is `∫ r/D = g + ∫ a/Ds` with `Ds` square-free. **The identity `g' + a/Ds = r/D` is checked exactly before use.**
3. **Logarithmic part.** Factor `Ds` into irreducibles `P_i`. Then `a/Ds = Σ A_i/P_i`, where `A_i` comes from `s·Q_i + t·P_i = a` and `Q_i = Ds/P_i`. Per factor:
   - `deg P = 1`: `k/(αx+β)` → `(k/α)·ln|αx+β|`
   - `deg P = 2`, numerator `Bx + C`: `B/(2a)·ln|P|` plus `K·∫1/P` with `K = C − B·b/(2a)` and `D = b² − 4ac`:
     - `D < 0`: `(2/√−D)·atan((2ax+b)/√−D)`
     - `D > 0`: `(1/√D)·ln|(2ax+b−√D)/(2ax+b+√D)|`
   - `deg P ≥ 3`: Lazard–Rioboo–Trager (`logpart.c`, `subres.c`); see §11.3.

```
∫ 1/((x-1)(x^2+1)) = 1/2 * ln(abs(x - 1)) - 1/4 * ln(abs(x^2 + 1)) - 1/2 * atan(x)
∫ 1/(x^2+1)^2      = 1/2x / (x^2 + 1) + 1/2 * atan(x)
∫ 1/(x^2-2)        = 1/4 * sqrt(2) * ln(abs((x - sqrt(2)) / (x + sqrt(2))))
```

The formulas here write `ln|u|`; the program prints it as `ln(abs(u))`, which nests and is easy for other tools to parse.

```
∫ 1/(x^3-2)        = RootSum(108t^3 - 1, t -> t * ln(x - 6t))
∫ x/(x^4+1)        = 1/2 * atan(x^2)
```

### 11.3 Lazard–Rioboo–Trager (`logpart.c`, `subres.c`)

For `∫ a/d` with `d` square-free, `deg a < deg d` (Bronstein, *Symbolic Integration I*, §2.5–2.8):

1. **Resultant.** `R(t) = res_x(d, a − t·d')`. Its roots are exactly the coefficients of the logarithms.
2. **Square-free split.** `R = Π Q_i^i`. Each nonconstant `Q_i` contributes `Σ_{Q_i(α)=0} α·log S_i(α, x)`. Here `S_i` is `d` when `i = deg d`, and otherwise the subresultant of degree `i` of `(d, a − t·d')`.
3. **Lazard–Rioboo normalization.** For the square-free parts `A_j` of `lc_x(S_i)`, divide `S_i` by `gcd(A_j, Q_i)^j`. This guarantees `S_i(α, x)` has the right degree at every root `α`.
4. **Per irreducible factor `M` of `Q_i`**, reduce `S_i` mod `M` and make it monic in x:
   - `deg M = 1`: `α·ln|S(α, x)|` with `α` rational.
   - `deg M = 2`, `α = p ± r√m` real: `p·ln|U² − mW²| + r√m·ln|(U + √m·W)/(U − √m·W)|`, with `S = S0 + t·S1`, `U = S0 + p·S1`, `W = r·S1`.
   - `deg M = 2`, `α = p ± r√m·i`: `p·ln|U² + mW²|` plus the arctangent part. For `m = 1` that part is `r·LogToAtan(U, W)`. Rioboo's recursion turns `i·log((A+iB)/(A−iB))` into a sum of `atan(poly)` terms that are continuous on the real line. Otherwise it is `−2r√m·atan(√m·W/U)`.
   - `deg M ≥ 3`: an exact `RootSum(M, t -> t * ln(S(t, x)))`.

**Subresultant chain.** `rt_chain` (`subres.c`) computes every `S_j` of `d` and `g = a − t·d'` over Q[t] in one subresultant pseudo-remainder sequence (Ducos 2000). With `A = S_d` regular (principal coefficient `s`) and `B = S_{d−1}` of degree `e`, δ = d − e: Lazard's shortcut gives `S_e = (lc(B)/s)^(δ−1)·B` (the `S_j` in between are zero), and `S_{e−1} = prem(A, −B) / (s^δ·lc(A))`. Every division is exact in Q[t] and asserted. The resultant is `S_0`. `test_subres_ref.c` keeps the definition (Sylvester-submatrix determinants evaluated at integer `t` and interpolated) as a test oracle, and `test_prop_subres_chain` checks the two agree exactly, including defective chains such as `1/(x^6+2)`. The chain is 6–90× faster than the determinants on degrees 8–25. For dense denominators, `reduce_monic`'s inverse of `lc(S)` modulo `M(t)` and the square-free gcds were the remaining cost. They are now modular (§9.3), for a 13–15× speedup at degree 30.

A final pass merges like terms (`c1·f + c2·f`), so cancelling arctangents disappear: `∫ x^2/(x^6+1) = 1/3 * atan(x^3)`.

**How it is verified.**
- When there is no `RootSum`, the test suite differentiates every antiderivative symbolically and compares it exactly with the integrand.
- A `RootSum`'s derivative at a rational point `x0` is `Σ α·S_x(α,x0)/S(α,x0)`. That equals the trace `Tr(t·S_x·S⁻¹ mod M)`, computed exactly from Newton power sums, so `RootSum` answers are also compared exactly at six rational points.
- `crosscheck_sympy.py` rechecks every `RootSum` answer at 60 digits.

### 11.4 Expanding RootSum when the roots are known (`rootform.c`, `surd.c`)

`RootSum` is exact but unfamiliar, so `integrate()` rewrites it into real `ln`/`atan` terms with radicals whenever the roots of the denominator factor `P` are known exactly. `integrate_rootsum()` skips this step and is what the exact tests use.

**When.** After LRT, only for a factor whose result contains a `RootSum`, so compact LRT answers such as `∫ x/(x^4+1) = 1/2 * atan(x^2)` are kept. `P` must be one of two shapes:
- **binomial** `p·x^n + q`: the roots are `ρ·e^(iπj/n)` with `ρ = |q/p|^(1/n)`, and `j` even or odd by the sign of `−q/p`;
- **rescaled cyclotomic in `x^s`**: `P(x) = Q(x^s)` with `Q = q_m·b^m·Φ_N(y/b)` and `b` rational. The roots are every `s`-th root of every `b·ζ_N^k` with `gcd(k, N) = 1`. Examples: `x^5 − 1`, `x^6 + 1`, `x^4 + 2x^2 + 4 = 4Φ_3(x^2/2)` (from `x^6 − 8`), and `x^6 + 3x^3 + 9`.

Anything else stays `RootSum`. For example, `x^5 + x + 3` has Galois group S5, so its roots have no radical form, and the Capelli factors `x^4 ± 2x^2 + 2` of `x^8 + 4` hold only some of the roots of a binomial.

**How.** For `h = a·(P')⁻¹ mod P`, `∫ a/P = Σ_r h(r)·ln(x − r)`:
- **Real root `r = ±ρ`:** gives `h(r)·ln|x − r|`.
- **Conjugate pair `r = ρe^(±iθ)`, `0 < θ < π`:** with `h(r) = H_re + i·H_im`, gives `H_re·ln(x² − 2ρcosθ·x + ρ²) − 2·H_im·atan(cscθ/ρ·x − cotθ)`. The atan argument is a polynomial in `x`, so it is continuous and needs no Rioboo step.
- **`H_re` and `H_im`:** both are `Σ_k h_k·ρ^k·cos/sin(kθ)`, sums of exact monomials.

**Exact constants (`surd.c`).** A `Mono` is `c · m^e · sqrt(s) · trig(angle·π)`:
- `R^e` is put in canonical form by extracting perfect powers, so `(1/1492992)^(1/6) = 2^(1/6)/12`.
- Trig values at multiples of π/6 and π/4 come from a table, so `cos(π/6) = sqrt(3)/2`. Other angles stay symbolic, folded into `(0, π/2)`, and print as `cos(2*pi/5)`. They are an `EXPR_TRIG` node. Values such as `cos(π/5) = (1 + √5)/4` are kept as `cos(pi/5)` on purpose, because it shows where the value comes from: the roots of `x^5 + 1` lie at angles kπ/5.
- **Square roots merge.** A square root goes into the radical slot when it is empty or already a square root, and squares move into `c`: `sqrt(2)·sqrt(3) = sqrt(6)`, and `sqrt(6)·sqrt(3) = 3·sqrt(2)` (for `x^6 + 216`). Beside a radical that is not a square root it stays separate, as `2^(1/12)·sqrt(3)`, unless it has the same base: `2^(1/4)·sqrt(2) = 2^(3/4)`.
- Like monomials are merged.

**How it is verified.**
- The C tests can't recognize radical and trig forms exactly, so they check these answers in double precision instead. The derivative is a 5-point stencil, taken only at points where two step sizes agree, which skips points near poles.
- `test_reg_integrate_radicals` pins five outputs.
- `test_prop_integrate_radicals` draws random `c·x^n ± e` and `x^(sN) ± b^N` denominators.
- `test_prop_integrate_derivative` checks `integrate_rootsum` exactly and `integrate` numerically.
- `crosscheck_sympy.py` adds a stream of binomial and cyclotomic integrands and checks them at 60 digits.

### 11.5 Definite integrals (`definite.c`, `constfold.c`)

`./polysolver integrate "f" a b` with rational bounds, or `inf`/`-inf` (also `oo`/`-oo`):

1. **Divergence.** If the denominator has a root in `[a, b]`, the integral diverges (every pole of a rational function has order ≥ 1). This is decided exactly: zero tests at `a` and `b`, then isolating intervals refined until each root is clearly inside or outside. An infinite side imposes no condition on the roots there. With an infinite bound, the integral converges only when `deg den ≥ deg num + 2`, since `f ~ c·x^(deg num − deg den)` far out.
2. **Antiderivative.** `F = integrate(f)`. If it contains a `RootSum`, only the numeric value is given: its complex logarithms have branch cuts, so `F(b) − F(a)` would not be trustworthy.
3. **Evaluation.** `expr_at` substitutes `x = a` and `x = b` and folds the constants exactly (`constfold.c`):
   - products collect their rational factor, and integer square roots multiply into one, with squares pulled out: `sqrt(2)·sqrt(6) = 2·sqrt(3)`;
   - sums collect their rational part and merge like terms;
   - `ln|c·m^e·…|` splits into logs of rationals, which merge when their coefficients are integer multiples (`-ln(3) + 2ln(2) = ln(4/3)`);
   - `atan` takes its exact values at `0`, `±1`, `±sqrt(3)`, `±1/sqrt(3)`, `±(sqrt(2) ∓ 1)` (π/8, 3π/8) and `±(2 ∓ sqrt(3))` (π/12, 5π/12), and `atan(−u) = −atan(u)`;
   - **Q(√m) arithmetic** (`qsurd.c`): a constant `p + q·sqrt(m)` with rational `p`, `q` is parsed into a `QSurd`, which supports `+ − · /`, an exact sign, and integer square roots in Z[√m]. In a folded sum, ln and atan terms whose arguments lie in one Q(√m) and whose coefficients share the same non-rational factor combine pairwise:
     - Coefficients in a ratio `p/q`, with `|p|, q ≤ 64`, merge as `(c/q)·(q·fn(u) + p·fn(v))`.
     - Logs: `q·ln|u| + p·ln|v| = ln|u^q·v^p|`. Surd logs merge only at integer ratios (`q = 1`). A fractional ratio would give arguments like `(1+√2)³ = 7+5√2`, whose cube root is not detected.
     - Atans add one at a time, tracked exactly as an angle `atan(w) + n·π/2`: `atan(w) + atan(v) = atan((w + v)/(1 − wv))`. This adds `sign(w)·π` when `1 − wv < 0`, and gives exactly `sign(w)·π/2` when `1 − wv = 0`. So `atan(1+√2) − atan(1−√2) = π/2` and `2·atan(1/3) = atan(3/4)`. The total angle fixes `w` and `n`, so a merged sum prints one way.
   - **Canonical log arguments.** Every `ln|g|` with `g` in Q(√m) is rewritten so that `g > 1` and has coprime integer parts, then reduced further:
     - a factor `√m` is split off while `m` divides `p`: `2+√2 = √2·(1+√2)`;
     - when `p` or `q` is negative, the conjugate is used instead: `ln g = ln|N| − ln|ḡ|`, with `N = p² − q²m`. So `ln(−1+3√2)` becomes `ln(17) − ln(1+3√2)`, and conjugate arguments print alike;
     - a square root is taken while `g` is a perfect square in Z[√m]: `3+2√2 = (1+√2)²`;
     - rational logs are written with `v > 1` and with perfect powers taken out: `ln(1/17) = −ln(17)`, `ln(485809) = 2·ln(697)`.

     Each rewrite splits off rational logs. Then terms are sorted: logs, atans, the rest, then multiples of π, ordered by printed form within each group. So `∫₀¹ 1/(x⁴+1) = 1/4·√2·ln(1+√2) + 1/8·√2·π`, and the integral over `[−1, 0]` prints the same.

     Check: for even integrands, `[a, b]` and `[−b, −a]` print identically in all 452 tested pairs.

     Surd logs are merged only with surd logs, and rational logs only with rational ones, since a merged surd log splits off rational logs and mixing the kinds would cycle. Merging ends because:
     - nothing creates an unnormalized log or a new surd log;
     - each surd merge removes a surd log;
     - each rational-log or atan merge removes a term of its own kind and adds no other.

     Pairs are tried in both orders, since a surd-log merge needs `c_j/c_i` to be an integer. The printer drops `abs` for a constant argument it can prove positive.

     Values in two radicals at once, such as `cos(π/12) = (√6+√2)/4`, or in a radical other than a square root, such as `2^(1/4)`, are left unmerged.
4. **atan jumps.** LRT can produce `atan(N/(c·W))`, for example `∫ (x²+1)/(x⁴+3x²+1) = 1/5·sqrt(5)·atan((x² − 1)/(sqrt(5)·x))`. This jumps by `±π` where `W` crosses zero. Each root of `W` inside `(a, b)` is isolated exactly, and the jump's direction is read from the sign of the argument on either side. The correction is `∓π` times the term's coefficient. A root of `W` at a bound uses the one-sided limit `±π/2`.
5. **Infinite bounds.** A term's limit comes from the leading behaviour of its argument, `u ~ L·x^d`. `expr_lead` in `constfold.c` finds it structurally: it handles `+ − · /` and integer powers, folds x-free subtrees as constants, and fails if the leading terms cancel. Then:
   - `c·atan(u)`: tends to `±c·π/2` when `d > 0`, with sign `sign(L)·(±1)^d`. It tends to `c·atan(L)` when `d = 0`, and to 0 when `d < 0`.
   - `c·ln|u|`: tends to `c·ln|L| + c·d·ln|x|`. The `ln|x|` parts sum to `(Σ cᵢdᵢ)·ln|x|`, which must vanish because the integral converges. So each log contributes `c·ln|L|`, and the logs of `1/(x²−1)` on `[2, ∞)` leave `1/2·ln(3)`.
   - A rational term from Hermite reduction is proper, so it tends to 0.
   - The atan jump count uses every root of `W` on the infinite side.
6. **Check.** The exact value is evaluated in double precision and compared with adaptive Simpson quadrature of `f`. The range is first cut at the real critical points of the denominator, the roots of `den'` isolated exactly. A complex pole close to the real axis makes a narrow peak at such a point. Without the cut, uniform pieces could step over the peak, or the half-line map could squeeze it to nothing: `1/((x−1000)²+1)` on the whole line was the failing case. Each piece is integrated in two halves, and each half has `f` Taylor-shifted exactly, in rationals, to its own end, so `f(o + u)` is evaluated with `u` small. Near the peak of `1/((x−50)² + 10⁻⁴)`, the expanded `x² − 100x + 2500.0001` would cancel catastrophically in double precision. Shifted, it becomes `u² + 1/10000`. A half-line beyond the outermost cut is mapped onto `t ∈ [0, 1]` by `u = ±(1−t)/t`, `du = dt/t²`. This is smooth because `|den|` is monotone there. It is bounded because `f = O(1/x²)`, and its value at `t = 0` is `lim x²f(x)`. With no cut, the whole line is split at 0. Every window starts from 64 pieces. The tolerance is relative 1e-10 of the summed rough estimate over all windows, not of each window, so a window where `f` is tiny is not driven into round-off. One budget of 4,000,000 evaluations is shared by all windows. Each window converts its numerator and denominator to doubles once (`DPoly`, Horner with the same rounding as `poly_evaluate`), not on every evaluation. The shift is not always the better basis: x^100 shifted to 1 is (1 + u)^100, whose binomial coefficients cancel catastrophically. So each polynomial is evaluated either shifted or in its own coefficients at o + u, whichever has the smaller Horner error bound `Σ|c_i|·|t|^i` at the window's far end (`Near`), and the shift is computed only up to degree 32 (it costs O(deg²) big rationals; x^99999 on [0, 1] used to hang). A finite window wider than 64 peak scales is cut geometrically from its origin, [0, s], [s, 2s], [2s, 4s], …, where s is half the least `|c_0/c_k|^(1/k)` of the denominator expanded at the origin, a lower bound on the distance to its nearest root: `1/(x²+1)` on [0, 1e40] and `1/(x²+1e-300)` on [0, 1] were lost among uniform pieces. When the exact value itself is beyond double range (`x²` on [0, 1e200]), there is nothing to compare, and the value is given with a note. If the exact and numeric values disagree by more than 1e-8 (relative), only the number is reported, with a note. The test suite asserts that this never happens on generated inputs.

   Roots at a finite end of the range are divided out before isolation (`roots_in`). Bisection mids are dyadic, so a root such as `−2/3` sitting exactly on a bound would never be separated from it.

Example: `∫_{-1}^{1} (x²+1)/(x⁴+3x²+1) dx = 1/5 * sqrt(5) * pi`. The naive `F(1) − F(−1)` would be 0, because the atan jump at `x = 0` would be missed.

Infinite examples: `∫_{-∞}^{∞} 1/(x⁴+1) = 1/2 * sqrt(2) * pi`, `∫_0^∞ x/(x³+1) = 2/9 * sqrt(3) * pi`, `∫_{-∞}^{∞} 1/(x⁶+2) = 1/3 * 2^(1/6) * pi`.

**Exponential sums and bounds q·π** (`definite_exp`). The same steps apply when f has exp, sin or cos, or a bound is a rational multiple of π (`pi`, `-pi/2`, `3pi/4`, parsed by `parse_endpoint`). f must be an exponential sum `c_0 + Σ c_w·e^w` (§11.7: no `ln`, no sum in a denominator). `expsum_shape` reports c_0, the product of the coefficient denominators, and, at each end, whether every `w ≠ 0` has `Re w → −∞`. The differences from the rational case:
1. **Poles.** A root of a coefficient denominator in `[a, b]` is not necessarily a pole, because the terms can cancel: `sin(x)/x` and `(x·cos(x) − sin(x))/x²` are finite at 0. Such a range is reported as unsupported. When f = c_0 is a rational function, the root is a real pole, and divergence is reported as before.
2. **Infinite ends.** An end at ±∞ is accepted only when every `w ≠ 0` has `Re w → −∞` there.
   - Then each term `y_w·e^w` of F tends to 0, since its coefficients are rational and `y_w ≠ 0` only for those w.
   - c_0 must decay like `1/x²`, as in step 1. Otherwise F's rational part is unbounded while the rest tends to 0, so the integral diverges.
   - An oscillating term with a decaying coefficient, such as `cos(x)/x` on `[1, ∞)`, gives a conditionally convergent integral. Quadrature cannot confirm it, so it is reported as unsupported. So is a growing term (`exp(x)` on `[0, ∞)`): no divergence is claimed there.
3. **q·π.**
   - In the span, q·π is replaced by `q·PI_50`, where `|π − PI_50| < 10⁻⁵⁰`. A polynomial with rational coefficients never vanishes at q·π, since π is transcendental.
   - So it is enough that the band `q·PI_50 ± (|q| + 1)·10⁻⁴⁹` holds no root of the pole polynomial, of an atan argument's denominator (for jumps), or of the other bound (for the order). `clear_of_bands` checks this exactly, and the result is reported as undecided otherwise.
4. **Evaluation.** `expr_at_pi` evaluates F at q·π exactly (`constfold.c`):
   - a polynomial p becomes `Σ c_k·q^k·π^k`, so π lands in the multiple-of-π slot and `π^k` becomes a factor;
   - `sin/cos(r·π)` takes `surd.c`'s table value (`trig_value`: `cos(π/3) = 1/2`, `sin(3π/4) = sqrt(2)/2`), and stays symbolic off the table (`sin(pi/5)`);
   - `sin/cos(a + r·π)` with 2r an integer is a quarter-turn shift: `sin(a + π) = −sin(a)`, `cos(a + π/2) = −sin(a)`;
   - other arguments stay symbolic (`cos(pi^2)`, `exp(pi)`);
   - at rational points, `exp(0) = 1`, `sin(0) = 0`, `cos(0) = 1`, `sin(−r) = −sin(r)` and `cos(−r) = cos(r)`.
5. **Antiderivative.** F comes from `integrate_exp`, or from the rational integrator when f = c_0. Without an elementary antiderivative, or when the integrator does not reach f, only the number is given: `∫₀¹ exp(−x²) ~ 0.7468241328`.
6. **Check.** Quadrature runs on `expr_eval(f)`: one window for a finite range, or half-lines mapped onto `[0, 1]` from the finite end (from 0 for the whole line). The value at `t = 0` is 0 because f decays exponentially. The exact value must agree to 10⁻⁸, relative to the larger of the value and the rough `∫|f|`. F(b) − F(a) can cancel to 0 while f reaches hundreds: D(Σ r·e^{−x}·sin/cos) on `[−π, ∞)` was the failing case.

Examples: `∫₀^π x·sin(x) = pi`, `∫₀^π x²·sin(x) = pi^2 − 4`, `∫₀¹ x·exp(x) = 1`, `∫₀^∞ x·exp(−x)·sin(x) = 1/2`, `∫₀^π cos(x + 1) = −2·sin(1)`, `∫₀^π 1/(x − 4) = ln|4 − π| − 2·ln(2)`.

Tests:
- `test_reg_definite_exp` pins 19 values, 2 numeric-only values and 6 refusals.
- `test_prop_definite_exp` integrates `D(Σ r(x)·exp(−x)·sin/cos(kx))` over random bounds `q·π` (q a multiple of 1/4) and `[q·π, ∞)`. It checks that each value equals `F(b) − F(a)` numerically and that every sin/cos is folded away.
- The SymPy `check_definite_exp` stream compares with `mpmath.quad` to 28 digits on rational, q·π and infinite bounds, and also checks the refusals at infinity. In 300 cases (seeds 100–399) there were 0 failures: 201 exact values, 42 numeric-only (with `exp(−x²)`), and 57 refusals at a non-decaying infinite end.

**A sum of exp or sin/cos in a denominator** (`definite_tower`, `deftower.c`). When `expsum` rejects a divisor with several terms, f is converted to `N/D` in one tower, θ = exp(v) or θ = tan(H) (§11.6, §11.8). Only finite ranges are taken: F may have infinitely many crossings on an infinite one, and terms such as `x − ln(exp(x) + 1)` cancel there. Then:
1. **Supported shape.** The coefficients of D must be constants, and the coefficients of N must have x-denominators with no root in `[a, b]`.
2. **tan(H).** H must be linear, `h1·x + h0`. All special points are of the form `α + β·π`, with `α = −h0/h1`.
   - **θ = ∞** at `H = π/2 + kπ`. f has a pole there when `deg N > deg D`. Otherwise these are the crossings where F may jump.
   - **θ = 0**, when `D(0) = 0`, at `H = kπ`. For k ≠ 0 the point is transcendental, so it is a pole (see the argument below). For k = 0 the point α is rational. There N is evaluated exactly: if it is nonzero, the point is a pole; if it is zero, the case is unsupported. `(x+1)/sin(x+1)` on `[−2, 0]` is finite.
   - These points are compared with the bounds exactly when their π parts agree, and otherwise with π to 50 digits. That leaves no tie at bounds such as π.
   - **Other real roots r of D** are at `H = atan(r) + kπ`. They are located in double precision from isolating intervals of width 10⁻¹⁵. A point within a relative 10⁻⁹ of a bound is reported as undecided.
3. **exp(v).** Only positive roots r of D matter, since `exp(v) > 0`.
   - r = 1 gives the rational point `−v0/v1`, tested exactly as above: `x/(exp(x) − 1)` on `[−1, 1]` is unsupported, not divergent.
   - Other roots are at `(ln r − v0)/v1`, located in doubles.
   - A nonlinear v with any positive root is unsupported.
4. **Why a root of D is a pole.** For algebraic r ∉ {0, 1}, `atan(r) + kπ` and `ln(r)` are transcendental by Lindemann–Weierstrass: `e^{2i·atan r} = (1 + ir)/(1 − ir)` and `e^{ln r} = r` are algebraic. So the point x is transcendental.
   - `N(·, r)` is a nonzero rational function with algebraic coefficients. If it vanished identically, the minimal polynomial of r over Q(x), which is its minimal polynomial over Q, would divide both N and D, but they are coprime.
   - So `N(x, r) ≠ 0`, and f has a pole of order at least 1, which is not integrable. The same holds for θ = 0 with k ≠ 0 (x is a nonzero rational multiple of π plus α), and for θ = ∞ with `deg N > deg D` (there `lc(N)` is a nonzero rational function at a transcendental x).
   - Quadrature cannot catch a wrong divergence verdict, so this argument is what such a verdict rests on.
5. **Jumps.**
   - Only `atan(P(tan H))` terms of F jump at a crossing. Where f is finite, the logs' `n·ln|θ|` parts cancel with `ln|cos H|`, so they do not jump.
   - Each crossing inside the range adds the jump read from the argument's sign on either side, as for the atan poles in step 4.
   - At a bound that is a crossing, atan terms take `±π/2` from the inside, and proper fractions in tan(H) tend to 0. A log of tan(H) or of cos(H) there leaves only the number.
   - Example: `∫₀^{2π} 1/(2 + cos x) = 2/3·√3·π`, because `2/√3·atan(tan(x/2)/√3)` drops by `2π/√3` at π. Over `[0, π]` the value is the one-sided limit, `√3/3·π`.
6. The quadrature check compares as before.

Tests:
- `test_reg_definite_tower` has 12 values, 8 poles (2 of them at a bound or with an offset argument), and 9 unsupported cases that must not be called divergent: `x/sin(x)`, `(x+1)/sin(x+1)` and `x/(exp(x)−1)` around their finite points, `(x+cos(x)−1)/sin(x)`, an infinite range, an x-dependent D, a nonlinear argument, a coefficient's pole, and a bound 10⁻¹³ below the pole of `1/(1 − tan x)` at π/4, which is undecided.
- `test_prop_definite_tan` checks `∫` of `1/(a + b·cos x + c·sin x)` over k full periods against `2πk/√(a² − b² − c²)`, and random `[q1·π, q2·π]` against quadrature.
- The SymPy `check_definite_tower` stream checks that 'diverges' is reported exactly when a numerically found real pole lies in the range, and checks the values against `mpmath.quad` to 28 digits. A double root such as `3 − 3·sin(x)` touches 0 without a sign change, so the checker confirms such minima with `findroot`. In 300 cases (seeds 100–399) there were 0 failures: 189 exact values, 35 numeric-only, and 55 poles, with 21 cases skipped where the denominator is near 0 at a bound.

### 11.6 Integration with exp and ln: the Risch algorithm (`expfield.c`, `risch.c`)

**Setting.** The field is Q(x)(θ), with θ = exp(v) and v a nonconstant polynomial, so θ' = v'·θ and θ is transcendental over Q(x). The integrands are Laurent polynomials `f = Σ a_k·θ^k` with `a_k ∈ Q(x)`.

`exp_form` builds the setting from the parse tree:
- It collects every `exp(u)` and canonicalizes each `u`, which must be a polynomial.
- It checks that all the `u` are rational multiples of the first argument w, and sets `v = w/l`, where l is the lcm of the ratios' denominators. Every `u` is then an integer multiple `k·v`, and v is oriented so that `lc(v) > 0`.
- It converts the tree into a fraction `f/den` over `K[θ]`, K = Q(x) (`kpoly.c`): `+ − · /`, integer powers capped at 256, and `exp(k·v) → θ^k`. The fraction is reduced by a gcd in K[θ] and `den` is made monic. For exp, θ is a unit of `K[θ, 1/θ]`, so θ powers move into the numerator and `θ ∤ den`. Sums and products of reduced fractions use Henrici's formulas:
- `a/b ± c/d = (a·d' ± c·b')/(b·d')`, with `g = gcd(b, d)`, `b' = b/g` and `d' = d/g`. It is reduced only by `gcd(numerator, g)`.
- A product cancels `gcd(a, d)` and `gcd(c, b)` before multiplying.
- Inversion needs no gcd.

So only small pieces ever meet a gcd. Re-reducing each result with one big Euclid over Q(x) has the same coefficient growth as over Q: on a degree-12 integrand it took 8.6 s, against 0.05 s now. When `den = 1`, f is a Laurent polynomial and the paths below apply directly. Otherwise the denominator steps further down run first.

Anything else is reported as unsupported, with a reason: independent exponentials, and `exp` of a constant or a non-polynomial.

**Theorem.** By Liouville's theorem, ∫ f is elementary iff `∫ f = ∫ a₀ + Σ_{k≠0} y_k·θ^k` with `y_k ∈ Q(x)` satisfying the Risch differential equation `y' + k·v'·y = a_k`. The solution is unique, because `y = C·θ^(−k)` is not rational. ∫ a₀ is the rational integration of §11.2–11.4.

**Solving the equation.** Here f = k·v' is a polynomial.
- **Poles.** Since `f` has no poles, at a pole p of y the term y' dominates, so `ν_p(y) = ν_p(a) + 1`. A simple pole of `a` therefore has no solution. This is the classic case: `∫ exp(x)/x` is not elementary.
- **Denominator.** Otherwise, with `D = den(a)`, the denominator of y is `h = gcd(D, D')`, which is `Π p^(e−1)` over `D = Π p^e`. `D | h²` is checked by an exact remainder.
- **Numerator.** `y = q/h`, where q is a polynomial with `h·q' + (f·h − h')·q = N·h²/D`. B = f·h − h' has degree `deg f + deg h`, which exceeds `deg(h·q') − deg q`. So the equation is triangular from the top: q is found term by term, and a nonzero remainder means there is no solution.

**Outcomes.**
- If every equation is solvable, the exp part is checked exactly: `d/dx Σ y_k θ^k = Σ_{k≠0} a_k θ^k`.
- If one equation is not, the result is a proof that no elementary antiderivative exists, and the failing term `a_k·exp(k·v)` is printed.
- Input outside the setting is "unsupported". It is never called non-elementary.

**Examples.**
- `∫ (x−1)·exp(x)/x² = exp(x)/x`: the denominator is h = x and q = 1.
- `∫ (2x²+1)·exp(x²) = x·exp(x²)`.
- `∫ exp(x²)`: here a = 1, h = 1 and B = 2x. The equation 2x·q + q' = 1 needs deg q = −1, so there is no solution.

**The ln (primitive) case.** Here θ = ln(u) with u ∈ Q(x) nonconstant, θ' = L = u'/u, and the integrand is a polynomial `Σ_{k=0}^{n} a_k θ^k`. Every `ln` argument must be `u^k` exactly for an integer k: `ln(x^2) = 2·θ` and `ln(1/x) = −θ`. The base is oriented so that `deg num ≥ deg den`. `ln(2x)` together with `ln(x)` would bring in the constant `ln 2`, so it is unsupported. Two independent logarithms, such as `ln(x)·ln(x+1)`, go to the two-level integrator (§11.9).

By Liouville's theorem, `∫ f = Σ_{k=1}^{n+1} b_k θ^k + (an integral over Q(x))`, with `b_{n+1}` constant and `b_k ∈ Q(x)`. Matching θ^k gives `a_k = b_k' + (k+1)·b_{k+1}·L`.
- **Top-down.** Write `b_k = r_k + d_k` with `d_k` constant. From k = n down to 1, `g_k = a_k − (k+1)·r_{k+1}·L` must satisfy `∫ g_k = r_k + (k+1)·d_{k+1}·θ`. This is *limited integration*: `hermite_split` gives `∫ g = R + ∫ A/D` with D square-free, and then `A/D = c·L` must hold exactly for a rational c. Any other log or atan at a level k ≥ 1 proves that ∫ f is not elementary.
- **Level 0** is free: `∫ (a_0 − r_1·L)` is an ordinary rational integral, with ln, atan or RootSum.
- **Check.** `d/dx Σ b_k θ^k + g_0 = f` exactly.

**ln examples.**
- `∫ ln(x)/x`: g = 1/x = 1·L, so `d_2 = 1/2` and the answer is `1/2·ln(x)^2`.
- `∫ ln(x)/(x+1)`: `A/D = 1/(x+1)` is not a multiple of `1/x`, so the integral is non-elementary (the dilogarithm).
- `∫ ln(x²+1)^2`: level 1 needs `atan(x)`, which is not a multiple of `ln(x²+1)`, so it is non-elementary.

**Denominators in θ** (`integrate_fraction` in `risch.c`). These cover `1/(1+exp(x))`, `1/(x·ln(x))` and `1/ln(x)`. For f = N/E, with E monic in θ and θ ∤ E for exp, the steps are:
- **Proper part.** For exp, write the numerator as `N·θ^(−m)`. Bezout gives `s·θ^m + t·E = N`, so `f = s/E + t·θ^(−m)`. In the exp case θ itself is *special* (`θ | θ'`), which is why it is split off first. The ln case has m = 0.
- **Hermite reduction over K[θ].** This is Bronstein's quadratic HermiteReduce:
  - The square-free parts V of E come from `d/dθ`.
  - Each step solves `B·U·D(V) + C·V = −A/j`, where D is the derivation: `D(θ) = v'θ` or `D(θ) = u'/u`.
  - Every factor of E is *normal*: gcd(V, D(V)) = 1, because the only special irreducible is θ for exp, and there are none for ln.
  - The result is `∫ b/E = g + ∫ h/d`, with d square-free and `deg h < deg d`.
- **Residue test (Rothstein–Trager).**
  - **The residues.** They are the roots of `res_θ(d, h − z·D(d))`, which up to a factor in K is the characteristic polynomial of `s = h·D(d)⁻¹ mod d`.
  - **The criterion.** ∫ f is elementary only if these roots are constants (Bronstein, *Symbolic Integration I*, Thm 5.6.1).
  - **Exact, without randomness.** That holds iff every power sum `Tr(s^k)`, `k = 1..deg d`, lies in Q. The traces use `Tr(θ^j)` from d's coefficients by Newton's identities.
  - **The proof.** A nonconstant power sum proves the integral is non-elementary, and the part `h/d` is printed.
  - **Otherwise** Newton's identities give `R(z) ∈ Q[z]`, which is factored over Q (`residue.c`).
  - **Rational roots.** For each rational root c, `g_c = gcd(d, h − c·D(d))`, scaled to a primitive polynomial `G_c` in `Z[x][θ]`, contributes `c·ln(G_c)`.
  - **Quadratic factors.** The roots are `α = a ± b·√m` (real) or `a ± b·√m·i` (complex). Write `M = m` or `M = −m`.
    - `g_α = gcd(d, h − α·D(d))` is taken over K(√M), where a polynomial is a pair `P + √M·Q` (`logpair_tmpl.h`). The x-denominators of P and Q are cleared together.
    - The pair `α·ln g_α + ᾱ·ln g_ᾱ` becomes `a·ln|P² − M·Q²|` plus:
      - for real α, `b·√m·(ln|P + √m·Q| − ln|P − √m·Q|)`;
      - for complex α, `−b·√m·LogToAtan(√m·Q, P)`.
  - **Rioboo with a scalar.** The complex case runs Rioboo's LogToAtan in K[θ] on states `2·atan(κ·A/B)`, with `κ = ρ·√m`:
    - if `B | A`, the term is `2·atan(κ·A/B)`;
    - if `deg A < deg B`, the state becomes `(−B, A)` with `κ' = 1/κ = √m/(ρ·m)`;
    - otherwise, with `B·D − A·C' = G`, the term is `2·atan(κ·(A·D + B·C'/κ²)/G)` and the state becomes `(D, C')` with the same κ.
    - So the arguments stay polynomial in θ, which means no jumps where P vanishes, and √m never enters the coefficient field.
    - The arctangents are checked exactly: `Σ 2ρ_i·D(T_i)/(1 + ρ_i²·m·T_i²) = 2(P·D(Q) − Q·D(P))/(P² + m·Q²)` in K(θ).
  - **Signs.** Both atan and the log ratio are odd, so arguments are written with a positive leading coefficient: `atan(exp(x))`, not `−atan(−exp(x))`.
  - **Degree 3 to 6** (`rootsum_tmpl.h`): for an irreducible factor f of degree n, the answer is `RootSum(f, t -> t·ln(g_t))`.
    - K(α) = K[t]/(f) holds polynomials in t over K of degree < n; f stays irreducible over K because its coefficients are constants. Products reduce mod f and the inverse comes from Bezout against f. `kpoly_tmpl.h` is instantiated over these elements for K(α)[θ].
    - The gcd runs against `e = gcd(d, Σ_j f_j·h^j·D(d)^(n−j) mod d)` (homogeneous Horner, all over K). It is the product of the n conjugate gcds, so each `g_α` has degree `deg e / n`. Both degrees are checked.
    - `g_α` is printed as `Σ t^i·L_i` with `L_i ∈ K[θ]` (the expression form of `EXPR_ROOTSUM`, `create_rootsum_terms`), monic in θ and uncleared like the pairs over K1.
    - **Exact derivative.** `Σ_α α·D(g_α)/g_α = Tr(α·D(g_α)·(e/g_α))/e`. `D` acts on each slice `L_i` (α is a constant), `e/g_α` is an exact division over K(α)[θ], and the trace is `Σ_i c_i·Tr(α^i)` coefficientwise, with the power sums of f's roots from Newton's identities. A wrong `g_α` leaves a denominator in the remainder, which is reported as an internal failure, never as an answer.
    - Degree 7 or more is unsupported (`exp(x)/(exp(7x)+2)`): the arithmetic over K(α) grows with n.
    - Example: `∫ exp(x)/(exp(3x)+2) = RootSum(108t³ − 1, t -> t·ln(exp(x) + 6t))`: the residues are `1/(3β²)` at the roots β of θ³ + 2, so `t³ = 1/108`, and `β = −6t`.
  - **Exact derivative.** Each term carries its derivative in K(θ). For a pair that is `a·D(N)/N + 2b·M·(P·D(Q) − Q·D(P))/(P² − M·Q²)`. These derivatives feed the remainder below.
- **Remainder.** `f − D(g) − D(log part)` is computed exactly and must be a polynomial (Laurent for exp) in θ. A failure here is an internal error. The remainder carries the exp correction `−Σ c·deg(G_c)·v'` and the `ln(content)` terms automatically. It is then integrated by the paths above, which check themselves.
- **Examples.**
  - `∫ 1/(1+exp(x)) = x − ln(exp(x) + 1)`: the residue is −1, and the remainder is `0 − (−1)·1·v' = 1`.
  - `∫ 1/(x·ln(x)) = ln(ln(x))`: s = 1.
  - `∫ 1/ln(x)`: here `s = 1/L = x` is not constant (this is the logarithmic integral).
  - `∫ x/(1+exp(x))`: the residue is −x (this is a dilogarithm).
  - `∫ exp(x)/(1+exp(2x)) = atan(exp(x))`: residues ±i/2, `g = θ + i`, and LogToAtan(1, θ) gives `2·atan(θ)` after one swap.
  - `∫ 1/(1+exp(x)+exp(2x)) = x − 1/2·ln(exp(2x) + exp(x) + 1) − 1/3·sqrt(3)·atan(1/3·sqrt(3)·(2·exp(x) + 1))`: residues `−1/2 ± i·√3/6`.
  - With θ = exp(x²), `x/(θ − 3)` has the constant residue 1/6 and *is* elementary, while `1/(θ − 3)` is not.

**Verification.**
- **Pinned cases** in `test_reg_integrate_exp`: 12 exact results, 5 proofs, 7 unsupported inputs, and parser cases.
- **Round trips** in `test_prop_integrate_exp`: integrating `d/dx F` returns F in canonical form. Each is also checked by a 5-point numeric derivative against the original input, which does not share the tower.
- **Non-elementary by construction:** `c·exp(v)/(x − r)` is added to a solvable integrand, including with r on an existing double pole, so that the remainder test decides.
- **SymPy:** the cross-check stream compares verdicts with SymPy's `risch_integrate`. For example, seed 20261113 gives 23 elementary results, which differentiate back exactly, and 16 proofs that SymPy confirms.
- **ln tests.**
  - `test_reg_integrate_log` pins 11 answers, 5 proofs and 5 unsupported inputs.
  - `test_prop_integrate_log` checks round trips with the numeric derivative, and that adding `c·ln(x)/(x − r)` is always rejected.
  - A separate SymPy stream for ln matched on all 400 integrands over 10 seeds: 257 elementary results and 143 proofs.
- **Denominator tests.**
  - `test_reg_integrate_tower_fractions` pins 22 answers (10 with irrational residues: atan, √3, and real √2 logs), 6 RootSums of degree 3 to 6 (each also checked by the numeric derivative), 6 proofs and 2 degree-7 cases.
  - `test_reg_integrate_tan` adds two RootSums in tan(x/2), of degree 4 and 6.
  - `test_prop_integrate_tower_fractions` integrates `(A/B^m + c·ln(G) + a·θ)'`, on 2 of 3 trials plus `D(U)/(U² + e)` (irrational residues), over five towers. B and G are linear or quadratic in θ, with coefficients in Q[x]. On 1 of 4 trials it also adds `D(θ)·(c1·θ + c0)/(θ³ + a·θ + b)`, whose residues are constants (11 of 60 trials give a RootSum). Each result is checked by the numeric derivative.
  - That test also adds a term with a nonconstant residue, which must be rejected.
  - The SymPy `check_frac` stream also adds `atan(U/√e)/√e` terms to half of its derivatives. It matched on all 400 integrands, seeds 100–499: 223 elementary results that differentiate back exactly, including the arctangent cases, and 177 proofs that SymPy confirms. None were unsupported.

### 11.7 sin and cos: exponential sums (`cfield.c`, `expsum.c`)

**Setting.**
- `sin u = (e^{iu} − e^{−iu})/2i` and `cos u = (e^{iu} + e^{−iu})/2`, for polynomial u.
- A product of exponentials is the exponential of the sum of the exponents.
- So every integrand built from `+ − ·`, integer powers, division by a single term, and exp, sin and cos of polynomials is a finite sum `Σ c_w·e^w`, with:
  - exponents `w = E + i·F` (E and F in Q[x]),
  - coefficients `c_w ∈ Q(i)(x)`, kept as `(re + i·im)/den` with a real monic `den`.
- `exp(x)·exp(x²)` is the single term `e^{x²+x}`. The same path takes exp-only integrands that are outside one tower, so `exp(x) + exp(x²)` is two terms.

**Theorem.**
- `D(c·e^w) = (c' + w'·c)·e^w` acts on each term on its own.
- Choose a basis of the Q-span of the exponents. The exponentials of the basis are algebraically independent (Risch structure theorem), and the terms are distinct monomials in them.
- So Liouville's theorem gives: ∫ is elementary iff `y' + w'·y = c_w` has a solution in Q(i)(x) for every w ≠ 0. The `w = 0` part is a rational integral.
- Exponents that differ by a constant c are not independent (`e^{i(x+1)} = e^i·e^{ix}`), but proofs stay exact:
  - e^c is transcendental (Lindemann–Weierstrass);
  - the image of `y ↦ y' + w'·y` is a subspace defined over Q(i);
  - so `c₁ + e^c·c₂` lies in it only if `c₁` and `c₂` both do. (Apply an automorphism of C that fixes Q(i) and moves e^c, then subtract.)
- Nonzero constant exponents (`exp(3)`, `cos(1)` alone) are reported as unsupported.

**Solving the equation over Q(i)** (`c_solve_rde`). This is the algorithm of §11.6 with complex leading coefficients:
- The denominator is `h = gcd(D, D')`, and `D | h²` is required (a simple pole has no solution).
- Then `h·q' + (w'·h − h')·q = N·h²/D` is solved term by term from the top. `lc(w'·h)` is a nonzero element of Q(i).
- Example: `∫ sin(x)/x` gives `c = ∓i/(2x)`, a simple pole, so there is no solution: it is the sine integral Si.

**Real output.**
- Terms are grouped as `E ± i·F`, with `lc(F) > 0`.
- The pair `c·e^{E+iF} + d·e^{E−iF}` is written as `e^E·(i(c − d)·sin F + (c + d)·cos F)`.
- Both coefficients must be real; a failure there is an internal error.
- Groups are ordered by F, then E. The rational part comes first: `∫ sin(x)^2 = 1/2x − 1/4 * sin(2x)`.

**Checks.**
- The answer is checked exactly: `D(Σ y_w·e^w) = f − c_0`.
- `simplify` prints the canonical form, which proves identities such as `sin(x)^2 + cos(x)^2 = 1`.

**Tests.**
- `test_reg_integrate_trig` pins 16 answers, 6 proofs, 6 unsupported inputs and 7 identities.
- `test_prop_integrate_trig` integrates `d/dx` of random `r(x)·exp(±x)·sin/cos(k·x or x²)` sums back to the same canonical form. The answer is unique because there is no `e^0` part. It also checks that adding `sin(x)/(x − 5)` is always rejected.
- The SymPy `check_trig` stream passes f rewritten in exponentials to `risch_integrate`. In 300 cases (seeds 100–399) there were 0 disagreements: 172 answers differentiate back to f exactly, 47 proofs of non-elementarity agree with SymPy, and in 81 more cases we prove non-elementarity where SymPy's Risch raises internally (these are unconfirmed).

### 11.8 sin and cos in a denominator: θ = tan(H) (`tanfield.c`, `risch.c`)

**Setting.** `expsum.c` divides only by single terms. When a divisor is a sum (`1/(1 + cos(x))`, `tan(x)`), `integrate_exp` switches to the hypertangent tower of Bronstein (*Symbolic Integration I*, §5.10):
- θ = tan(H) for a nonconstant polynomial H, with `θ' = H'·(1 + θ²)`. The base field stays K = Q(x), so `kpoly.c` and `residue.c` work unchanged.
- Every sin/cos argument must be an integer multiple k·B of one polynomial B.
- `cos(jH) + i·sin(jH) = (1 + iθ)^j·cos(H)^j`, and `cos(H)² = 1/(1 + θ²)`. So an odd j leaves one factor cos(H). `tanfield.c` computes with pairs `a + cos(H)·b`, a and b in K(θ). The integrand lies in K(θ) exactly when b = 0.
- **Choice of H (Bioche's rule).**
  - A parity pass gives sin(kB) and cos(kB) the parity k mod 2. Sums must agree, and products add parities.
  - Even parity means the integrand is unchanged by B → B + π, so H = B (`tan(x)` stays `tan(x)`).
  - Otherwise H = B/2. Then every j is even, so b = 0 always.
  - `Tower.half` records the choice.

**The special factor.** `1 + θ²` divides its own derivative. So the denominator is split as `D = (1 + θ²)^m·E` with `gcd(E, 1 + θ²) = 1`. Then, by Bezout, `f = P + S/(1 + θ²)^m + A/E`:
1. **Normal part A/E.** Hermite reduction and the Rothstein–Trager residue test, exactly as in §11.6.
2. **Remainder.** `r = f − D(g) − D(logs)` must lie in `K[θ, 1/(1 + θ²)]`, which is checked exactly. Its polynomial part is reduced from the top by `D(c·θ^(n−1)) = c'·θ^(n−1) + (n−1)·c·H'·(θ^(n−2) + θ^n)`, down to `p1·θ + p0`.
3. **Polynomial part.**
   - `p1·θ = c·D(ln(1 + θ²))` with `c = p1/(2H')`.
   - If c is not a constant, f is not elementary.
   - When it is, the answer gets `c·ln(1 + θ²) = −2c·ln|cos H|`.
4. **Special part.**
   - `p0 + S/(1 + θ²)^m`, with `θ^a/(1 + θ²)^m = sin(H)^a·cos(H)^(2m−a)`.
   - This is a trig polynomial, integrated by `expsum.c`: one RDE per exponent, and `∫ p0` by the rational integrator.

**Why the parts are independent.** Put `z = e^{2iH}`, so `θ = −i(z − 1)/(z + 1)`, an exponential tower over Q(i)(x):
- The special part is a Laurent polynomial in z.
- The polynomial part has its poles at z = −1. Its residue there is `−p1/H'`, which is the test in step 3.
- A proper fraction in z (z ∤ den) differentiates to a proper fraction.
- So only the special part contributes nonzero powers of z, and §11.7's theorem applies to them. Each non-elementary verdict is a theorem about one part.

**Output.**
- θ prints as `tan(H)`: the tree is `sin(H)/cos(H)`, and the printer shows it as tan.
- Logs of polynomials G(θ) of degree n can be rewritten: `ln|G(tan H)| = ln|Ĝ| − n·ln|cos H|`, where `Ĝ = cos(H)^n·G(tan H)` is a trig polynomial.
  - The rewrite is used only when all the ln|cos H| coefficients cancel. Then Ĝ is put in canonical form by `simplify_expsum`.
  - When H = B/2 and n is odd, Ĝ is squared first, with half the coefficient, so that it is written in B.
  - Example: `∫ cos(x)/(1 + sin(x)) = 2·ln|cos(x/2)| + 2·ln|tan(x/2) + 1| = ln|1 + sin(x)|`.
  - When the coefficients do not cancel, the tan form stays: `∫ 1/sin(x) = ln|tan(x/2)|`.
- This rewrite runs after the exact checks, so it is covered by numeric derivative checks of the printed answer (tests and SymPy).

**differentiate** uses the same tower only when H = B: d/dx tan(x) = `tan(x)^2 + 1`. In tan(x/2) the exact derivative of `1/(1 + cos(x))` is `1/2·tan(x/2)^3 + 1/2·tan(x/2)`, which is correct but unreadable, so that case is reported as unsupported.

**Tests.**
- `test_reg_integrate_tan` pins 19 answers, each also checked numerically in its printed form. It also has 5 proofs:
  - `x·tan(x)`, `x·sin(x)/(1 + cos(x))` and `x²/(1 + cos(x))` by p1;
  - `1/(x + sin(x))` by residues;
  - `sin(x)/x + tan(x)` by the special part.
- It also has 3 unsupported inputs and the derivatives.
- `test_prop_integrate_tan` builds `D(N/E) + c·D(G)/G + D(P)/(1 + P²)` from random trig polynomials, using the quotient rule by hand, so the input is in sin/cos of x. It checks the printed answer's derivative numerically, and that adding `x·sin(x)/(1 + cos(x))` is always rejected.
- The SymPy `check_tan` stream compares with `risch_integrate` on f in exponential form. SymPy is given 20 CPU seconds per case, and answers are checked numerically to 35 digits. In 300 cases (seeds 100–399) there were 0 disagreements: 204 answers differentiate back to f, and 96 proofs of non-elementarity agree with SymPy (none unconfirmed, none unsupported).

### 11.9 Two levels: ln over an exp or ln tower (`risch2.c`)

**Setting.** θ2 = ln(u) over K1 = Q(x)(θ1), where θ1 = exp(v) or ln(u1) as in §11.6, u ∈ K1, and f is a polynomial in θ2 with coefficients in K1.
- Examples: `exp(x)·ln(x)` (K1 = Q(x)(eˣ), u = x), `ln(ln(x))/x` (K1 = Q(x)(ln x), u = ln x) and `ln(x)·ln(x+1)` (K1 = Q(x)(ln x), u = x + 1).
- The path is entered only when the one-tower conversion fails on an integrand with ln. So no earlier answer changes. It serves `integrate` only: `differentiate` and `simplify` still refuse two towers.
- **Conversion.**
  - The outer argument is chosen among the ln arguments: nested ones (containing exp or ln) first, then the others from the last.
  - K1 is the tower of the rest of f, with each outer ln replaced by x (argument collection is syntactic), together with u. This matters because in `ln(ln(x))/x` the tower appears only inside u.
  - f is then converted coefficientwise in that fixed tower (`exp_form_in`).

**Levels** (Bronstein, *Symbolic Integration I*, §5.8; Liouville's theorem as for one ln, §11.6).
- `∫ f = Σ_{k=1}^{n+1} b_k·θ2^k + ∫ g_0`, with b_{n+1} constant, b_k ∈ K1, and `a_k = D(b_k) + (k+1)·b_{k+1}·η`, where η = D(u)/u.
- From the top, g_k = a_k − (k+1)·r_{k+1}·η must equal D(r_k) + c·η with r_k ∈ K1 and c ∈ Q. This is limited integration over K1 (Bronstein §7.2), and then b_{k+1} = r_{k+1} + c/(k+1).
- Level 0, `g_0 = a_0 − b_1·η`, is any integral over K1 and goes to the one-tower integrator.

**Limited integration over K1, by Hermite reduction alone.** Write g = D(R) + h/d + P, with d square-free and P a polynomial in θ1 (Laurent for exp).
- R′ ∈ K1 with only simple normal poles in D(R′) has no normal poles at all, so it is a polynomial. Hence the simple part of an element of D(K1) + Q·η is c times η's.
- η's polynomial part p_η lies in Q(x):
  - for exp, `D(G)/G` has only a θ1⁰ Laurent part;
  - for ln, it is proper in θ1.
- So g = D(r) + c·η exactly when:
  1. `h_g/d_g = c·h_η/d_η`, which fixes c when η has a simple part;
  2. the θ1^k parts of P_g with k ≠ 0 integrate in K1: for exp, an RDE per k (`exp_part`); for ln, the levels of §11.6 (`solve_log_levels`), which leave the θ1⁰ remainder q = p_0 − b_1·L1;
  3. at θ1⁰, `q − c·p_η` has an integral in Q(x), up to `d·L1 = D(d·θ1)` for a ln base. In Hermite simple parts over Q(x) this reads `s(q) = c·s(p_η) + d·L1`, solved for constants over Q (`linsys_solve`).
- Then `r = R_g − c·R_η + Y + d·θ1 + ∫(q − c·p_η − d·L1)`, and `g − D(r) − c·η = 0` is checked exactly.
- The residue test is not needed here. Cubic residues therefore cannot make a level "unsupported": either a level fails, which is a proof, or it is solved exactly.

**Transcendence.** θ2 is a new transcendental only if η ∉ D(K1), meaning not both h_η = 0 and s(p_η) ∈ Q·L1. The case `ln(x)·ln(2x)` (ln 2 is new) and `exp(x)·ln(exp(x))` (= x) is reported as unsupported.

**Verdicts.** For a polynomial in θ2, non-elementary is claimed in exactly three cases:
- a level's c-free part fails (step 2);
- no (c, d) exists (steps 1 and 3);
- level 0 is non-elementary over K1 (its own proof, §11.6).

Denominators in θ2 add two more cases, given below. Everything else is unsupported: shapes outside this setting, internal check failures, θ2 in K1, and residues of degree 7 or more. The final answer is checked exactly: `D(Σ b_k θ2^k) + g_0 = f`.

Examples:
- `∫ exp(x)·(ln x + 1/x) = exp(x)·ln|x|`, with c = 0.
- `∫ ln(ln x)/(x·ln x) = 1/2·ln|ln|x||²`, with c from the simple part.
- `∫ (eˣ + 1/x)·ln x + eˣ/x = 1/2·ln|x|² + eˣ·ln|x|`, with c from the θ1⁰ part.
- `∫ (1/x + 1/(x+1))·ln(x+1) + ln(x)/(x+1) = ln|x+1|·ln|x| + 1/2·ln|x+1|²`, which uses the d·L1 freedom.
- `∫ exp(x)·ln(x)` is non-elementary: level 0 is `−eˣ/x`, the exponential integral. So is `∫ ln(x)·ln(x+1)`, which leads to the dilogarithm.

**Denominators in θ2** (`integrate_frac2`). Every factor of a denominator is normal, because θ2 is a primitive monomial with η ∉ D(K1), and that is checked first. So the method of §11.6 carries over to K1[θ2]:
- Bezout splits off the polynomial part.
- Bronstein's quadratic Hermite reduction runs with the derivation `D(Σ b_k θ2^k) = Σ (D(b_k)·θ2^k + k·b_k·η·θ2^(k−1))`.
- The residues of the square-free part h/d must be constants. Their power sums `Tr(s^k)`, with `s = h·D(d)⁻¹ mod d`, are computed from Newton's power sums of d's roots over K1, and each must lie in Q.
- `R(z)` then follows over Q by Newton's identities (`from_power_sums`, shared with `residue.c`). It is factored with `factor.c`, and each rational root c gives `c·ln(gcd(d, h − c·D(d)))`.
- A quadratic factor gives the conjugate pair exactly as over Q(x): `g_α = P + √M·Q` over K1(√M)[θ2], then `a·ln|P² − M·Q²|` plus real logs `b·√m·ln|(P + √m·Q)/(P − √m·Q)|` or Rioboo's arctangents of polynomials in θ2, whose derivative is checked exactly. g_α stays monic in θ2, so P and Q may have coefficients that are fractions in θ1 (`atan(ln|x+1|/ln|x|)`); clearing them would put `2a·D(c)/c` into the remainder for no gain. Arctangents of a constant are dropped (derivative 0). Residues of degree 3 to 6 give a RootSum through `rootsum_tmpl.h` over K1, exactly as over Q(x); `g_α` may have coefficients in K1, as in `RootSum(108t³ − 1, t -> t·ln(ln|x+1| − 6·ln|x|·t))`.
- `f − D(g) − D(logs)` must be a polynomial in θ2. This is checked exactly, and the polynomial goes to the levels above.
- Verdicts:
  - Non-elementary is claimed when a power sum is not in Q (Bronstein Thm 5.6.1 for a primitive monomial), or when the polynomial remainder's levels fail.
  - Residues of degree 7 or more are reported as unsupported, as are denominators of degree above 6 and internal failures. Cubic ones, such as those of `1/(x·ln(x)·(ln(ln(x))³ − 2))`, give a RootSum.
- Examples:
  - `∫ 1/(x·ln x·ln ln x) = ln|ln|ln x||`, with residue 1;
  - `∫ 1/(x·ln x·(ln ln x + 1)²) = −1/(ln|ln|x|| + 1)`, by Hermite;
  - `∫ (1/x + eˣ)/(ln x + eˣ) = ln|ln|x| + eˣ|`;
  - `∫ 1/(x·ln x·((ln ln x)² + 1)) = atan(ln|ln|x||)`, with residues ±i/2, and `∫ ln ln x/(x·ln x·((ln ln x)⁴ + 1)) = 1/2·atan(ln|ln|x||²)`, where LogToAtan takes more than one step;
  - `1/(ln x·ln(x + 1))` is non-elementary with either ln on top.

**The pair layer.** `logpair_tmpl.h` holds the algebra of a conjugate pair once: the gcd over K(√M), `P² − M·Q²`, the pair's derivative, Rioboo's LogToAtan steps and the exact check of the arctangents. `residue.c` instantiates it over Q(x) (its output is byte-identical to the earlier `ksurd.c` version on the whole suite) and `risch2.c` over K1; only the printing differs. `rootsum_tmpl.h` uses the same macros (plus `LT_C_COPY`, `LT_C_IS_ZERO`) for residues of degree 3 to 6, and `logterm_undef.h` undefines them after both.

**The polynomial layer.** `kpoly_tmpl.h` holds divmod, gcd, Bezout, Yun's square-free decomposition, monic and d/dθ. They are written once over macros for the coefficient operations.
- `kpoly.c` instantiates them over QFn (K = Q(x)), and the output is unchanged, which the whole suite checks.
- `k2poly.c` instantiates them over KFrac (K1) for θ2. `test_prop_k2poly` checks the identities there: a = q·b + r, gcd·cofactor = input, s·a + t·b = 1, and that Yun's parts multiply back.
- Euclid over K1 reduces fractions in θ1 at every step, and its coefficients swell, which is why the degree cap exists. The cap is applied while f is converted, too: a sum or product whose denominator in θ2 would exceed it is refused before any gcd (a sum of four fractions of degree 2 used to take minutes there).
- `k2frac_add` uses Henrici's addition: with `g = gcd(d_a, d_b)`, the numerator `t = n_a·(d_b/g) ± n_b·(d_a/g)` can share factors only with g, so only `gcd(t, g)` is taken, which is trivial when the denominators are coprime.
- **Speed.** The template's `divmod` works in place on a coefficient window, and `gcd` and `inverse` keep every remainder monic, so the dividing step needs no inverse. QFn uses Henrici's formulas too: `qfn_add` takes `gcd(d_a, d_b)` and then only `gcd(t, g)`, and `qfn_mul` cancels `gcd(n_a, d_b)` and `gcd(n_b, d_a)` before multiplying. Both need no gcd at all when the operands are polynomials (denominator 1). Together with §5's list-free core this took the suite from 15.7 s to 1.9 s of CPU time, with byte-identical CLI output on the golden corpus.
- The pair gcd over K1(√M) runs against `gcd(d, (h − a·D(d))² − M·b²·D(d)²)`, the product `g_α·g_ᾱ` found over K1 without surds, rather than against all of d.

**Tests.**
- `test_reg_integrate_tower2` pins 25 answers, each also checked numerically (three with a Hermite part or a negative power in θ1, fifteen with θ2 in a denominator: two residues, a non-monic denominator, a zero residue after Hermite, and eight with quadratic residues: real and complex pairs, a ≠ 0, Q with a coefficient in K1, a two-step LogToAtan, and pairs beside rational residues), 4 RootSums (two with a coefficient in K1), 15 proofs (both orientations of `1/(ln(x)·ln(x+1))`, nonconstant residues over quadratics) and 5 unsupported inputs (one with residues of degree 7).
- `test_prop_integrate_tower2` integrates the derivative of random `Σ c_k(x)·{θ1, θ1², θ1/(θ1 + 2)}·θ2^k`, sometimes with `c·θ2²/(θ2 + 3)`, and sometimes plus `D(atan(r·θ2))` or `D(ln(θ2² + r))` (quadratic residues), and on 1 of 4 trials without another θ2 denominator `D(θ2)·(c1·θ2 + c0)/(θ2³ + a·θ2 + b)` (a RootSum over K1), in four tower pairs. It adds a non-elementary integrand of the same towers, and checks that the sum is always rejected.
- The SymPy `check_tower2` stream gives these towers to `risch_integrate` unchanged. Coefficients include `1/(θ1 + 2)`; 40% of the integrands have a term `c/(θ2 + c')`, 30% a term `c/(θ2² + c')`, and a fifth of the derivatives an atan or ln of a quadratic in θ2. In 300 cases (seeds 0–299) there were 0 disagreements: 155 answers check numerically, 134 proofs of non-elementarity agree with SymPy, 3 proofs remain unconfirmed because SymPy's Risch raised or timed out, and 8 are refused by the degree cap.

---

### 11.10 One square root (`algsqrt.c`)

**The field.** An expression with square roots of one polynomial U is an element `a + b·y` of Q(x)(y), y² = U, with a, b ∈ Q(x):
- Products use y² = U, and inverses use the conjugate: `1/(a + b·y) = (a − b·y)/(a² − b²·U)`. The norm is zero only for 0, because U is not a square.
- `sqrt(g)` normalizes g = c·w (w primitive, positive leading coefficient) to `r·sqrt(±m·w)` with |c| = r²·m, so `sqrt(4x+4) = 2·sqrt(x+1)`.
- Every root in an expression must give the same radicand.
- A constant square (`sqrt(4/9)`) is rational. A negative constant, or a perfect square (`sqrt(x^2) = |x|`), is refused.
- Half-integer powers `g^(n/2)` are `sqrt(g)^n`.

The derivative is `D(a + b·y) = a' + (b' + b·U'/(2U))·y`.

**Integration**, for deg U ≤ 2. `∫ a` goes to the rational integrator. For `∫ b·y`, the conic y² = U gets a rational parametrization x = X(t), y = Y(t), with its inverse t = T(x, y):

| U | t | X(t), Y(t) |
|---|---|---|
| b₁x + c | y | X = (t² − c)/b₁, Y = t |
| s²x² + bx + c | y + s·x (Euler) | X = (t² − c)/(2st + b), Y = (st² + bt + sc)/(2st + b) |
| a rational point (x₀, y₀) | (y − y₀)/(x − x₀) | X = (x₀t² − 2y₀t + u₁ − a·x₀)/(t² − a), Y = (−y₀t² + u₁t − a·y₀)/(t² − a), u₁ = U'(x₀) |

- **The rational point.** The search tries x₀ = p/q by increasing |p| + q, up to 40, with x₀ = 0 first, and needs U(x₀) to be a rational square.
- **No point.** A conic with no rational point, such as x² + y² = 3, has no parametrization over Q. It goes to the reduction below.
- **The exact check.** `X(T) = x` and `Y(T) = y` must hold in Q(x)(y). Then, for `F'(t) = R(t) = b(X)·Y·X'`, the chain rule gives `d/dx F(T) = R(T)·T' = b(x)·y`, because `X'(T)·T' = 1`.
- **Back-substitution.** R is integrated by `integrate()`. The answer's terms that are rational in t are summed in Q(t) and substituted once, so they come back as one normal form `a + b·y`, written over one denominator, with an additive constant dropped. The ln, atan and RootSum terms are substituted node by node; a RootSum in t becomes the expression form `Σ αⁱ·cᵢ(x, y)`.

**Examples.**
- `∫ 1/sqrt(x²+1)`: here t = x + y, X = (t² − 1)/(2t), Y = (t² + 1)/(2t), and R(t) = 1/t, so the answer is `ln|x + sqrt(x² + 1)|`.
- `∫ sqrt(1−x²)`: the point (0, 1) gives t = (y − 1)/x, and the answer is `x/2·y − atan((y − 1)/x)`. The argument tends to 0 as x → 0, so there is no jump.

**Conics without a rational point.** Write `b·y = R/y` with R = b·U and take R's partial fractions over Q. U is irreducible here (a rational root would be a rational point), so every other factor d is coprime to U.
- **Powers of U** (and the polynomial part) form `N/(U^m·y)`, reduced as below.
- **Other poles d^j, j ≥ 2: Hermite steps.** `d/dx(Q·y/d^(j−1)) = (Q'U·d + Q·U'·d/2 − (j−1)·Q·d'·U)/(d^j·y)`, so `Q ≡ −p/((j−1)·d'·U) mod d` lowers the power by one. The quotient of the last step goes to the polynomial part.
- **Simple poles r/(d·y)** give the logarithmic part (`sqrtlog.c`). Each root ρ of d contributes `c·∫ dx/((x−ρ)·y)` with `c = r(ρ)/d'(ρ)`, and for any w with w² = U(ρ)
  `∫ dx/((x − ρ)·y) = −(1/w)·ln((2U(ρ) + U'(ρ)(x − ρ) + 2w·y)/(x − ρ))`.
  With t = −c/(2w), the term is 2t·ln G, and s = t² = c²/(4U(ρ)) lies in Q(ρ):
  - **s rational** (always for linear d; also `1/((x²+1)·sqrt(3−x²))`): w = c/(2√s) at every root, so the terms combine through the norm `E = Π_ρ (L(ρ) + B(ρ)·ν)`, ν = y/√s, L = 2U(ρ) + U'(ρ)(x − ρ). E is computed by interpolation: determinants of multiplication matrices in Q[ρ]/d on an (n+1)² grid of (x, ν). With ν² = U/s, `E = P0 + P1·ν`, and the answer is `−2√s·(ln|P0 + P1·y/√s| − ln|d|)` for s > 0, or `2√|s|·atan(√|s|·P0/(P1·y))` for s < 0.
  - **s generating Q(ρ)**: R = minpoly(s), found from the powers of s in Q[ρ]/d by linear algebra, which also gives ρ, U(ρ) and c as polynomials in s. The answer is `Σ_t t·ln(L(ρ) + 2w·y)` over the roots t of R(t²), one RootSum per irreducible factor. The ln(x − ρ) parts cancel between t and −t.
  - **s generating a proper subfield** (`1/((x⁴+1)·sqrt(3−x²))`): unsupported.
- **Steps for m ≥ 1.** Use `d/dx(Q·y/U^m) = (Q'U + (½ − m)·Q·U')/(U^m·y)` and choose `Q ≡ N·U'^{-1}/(½ − m) mod U`. U is square-free, so U′ is invertible mod U. Then U divides the difference, and m drops by one.
- **Base case m = 0.** Write `N = Q'U + Q·U'/2 + λ`. Q is solved from its top coefficient down: x^j in Q puts `a·(j + 1)` on x^(j+1).
- **The last integral** is `∫ 1/y`:
  - `(1/√A)·atan(√A·(x + b/2a)/y)` with A = −a > 0, continuous inside the domain;
  - `(1/√a)·ln|2ax + b + 2√a·y|` when a > 0.
- **The exact check.** The algebraic part is verified in Q(x)(y): `D(G·y) + λ/y + Σ r_i/(d_i·y) = b·y`, G holding every Hermite step. The ln/atan/RootSum terms of the simple poles have surd coefficients outside Q(x)(y), so the whole answer is also checked by central differences at sample points inside the domain; a mismatch is reported as an internal failure, never returned.
- **Examples.**
  - `∫ 1/sqrt(3 − x²) = atan(x/sqrt(3 − x²))`.
  - `∫ x²·sqrt(5 − 2x²)` comes out as `(1/4x³ − 5/16x)·sqrt(5 − 2x²) + 25/32·√2·atan(√2·x/sqrt(5 − 2x²))`.
  - `∫ 1/((x − 2)·sqrt(3 − x²)) = atan((3 − 2x)/sqrt(3 − x²))` (U(2) < 0), and `∫ 1/(x·sqrt(3 − x²))` is a ln with √3.
  - `∫ 1/((x² + x + 1)·sqrt(3 − x²))` is a RootSum over t⁴ + 7/156·t² + 1/1872.

**Definite integrals** (`sqrtdef.c`), for rational and infinite bounds.
- **Exact range checks.**
  - U ≥ 0 on [a, b]: no root of U strictly inside (from refined isolating intervals), and U > 0 at the middle.
  - No pole of a or b strictly inside.
  - At each end e, the exact order of f = a + b·y must exceed −1. It is the smaller of the orders of a and b·y (half-integers at a root of U), unless their leading terms cancel (α = −β·sqrt(U(e)), checked exactly as α² = β²·U(e) with opposite signs). Then it is `ord(a² − b²U) − ord(a)`.
- **Evaluating F(b) − F(a)**, term by term.
  - Constant factors stay as they are.
  - ln and atan take their arguments' values.
  - Algebraic pieces, and those arguments, go through Q(x)(y) (`sqrt_value_at`).
  - **0/0 at an end.** Where the common denominator D vanishes at the end, `(A + B·y)/D = (A² − B²U)/(D·(A − B·y))`. The factor (x − e)^k is cancelled from `A² − B²U`, which handles cases like `(sqrt(1 − x²) − 1)/x` at 0.
  - **Mixed arguments.** An argument that is not in Q(x)(y) as a whole (it mixes `sqrt(3)` with y) is evaluated side by side through + − · /, each side through Q(x)(y) with its limit, and folded by `expr_at`, merging `sqrt(3)·sqrt(5)`. A quotient whose sides are both infinite is rescaled by (x − e)^k, k ≤ 4, distributed over the sums. Otherwise x is substituted everywhere, inside roots too.
  - **Algebraic terms together.** All summands without ln or atan are evaluated as one element of Q(x)(y), since each may be infinite at an end while the sum has a limit.
  - **atan jumps.** Where an atan argument `c·N/(W·y)` has a pole inside (a zero of W, isolated exactly), the jump ±π is subtracted, with the signs read just beside it.
  - **Infinite atan arguments.** An atan whose argument tends to ±∞ at an end (`atan(y/x)` at 0) takes ±π/2. The sign comes from the argument just inside the range, which must grow toward the end.
  - The constants are summed by `const_fold_sum`.
- **The numeric check.** The value is compared with tanh-sinh quadrature, done on each half from its own end. The polynomials of `a + b·sqrt(U)` are shifted exactly there, so u near 0, a singular end, keeps its precision: sampling `1 − d` near 1 would round away a √d-sized sliver. Without an exact value, the number is given with the reason.
- **The conjugate in quadrature.** Where `a` and `b·y` nearly cancel, the integrand is evaluated as `(a² − b²U)/(a − b·y)`, with the exact powers of u divided out of every shifted polynomial so that u^12/u^12 does not underflow.
- **Infinite bounds.** x = 1/u maps [c, ∞) to (0, 1/c], with `sqrt(U(1/u)) = sqrt(V(u))/u^k`, V = u^(2k)·U(1/u), k = ⌈deg U/2⌉. A −∞ end is mirrored by x → −x first, and (−∞, ∞) is split at 0.
  - Convergence is the exact order test at u = 0 for G(u) = f(1/u)/u², so `(sqrt(x² + 1) − x)²` converges and `sqrt(x² + 1) − x` is refused.
  - The number is the quadrature of G on [0, 1/c] (plus [c, 1] when c ≤ 0).
  - The exact value is `lim F − F(c)` for one antiderivative F: the limit is F(1/u) at u → 0⁺ with the evaluation above, minus the atan jumps on (c, ∞) (roots isolated up to Cauchy's bound). Before the limit, ln terms with opposite equal coefficients are merged into one ln of a quotient, and products are split into their x-free factors and one element of Q(x)(y), so that ∞·0 inside a rescaled quotient still has its limit.
- **Results.**
  - `∫₀^∞ 1/(1 + x²)^(3/2) = 1`, `∫₀^∞ (sqrt(x² + 1) − x)² = 2/3`, `∫₀^∞ 1/((x + 1)·sqrt(x)) = π`.
  - `∫₋₁^½ 1/((x² + 1)·sqrt(3 − x²))` crosses the atan jump at 0.
  - `∫₀¹ sqrt(1 − x²) = π/4`, `∫₀¹ 1/sqrt(x(1 − x)) = π` and `∫₀¹ 1/sqrt(1 + x²) = ln(1 + √2)`.
  - `∫₀¹ 1/sqrt(x³ + 1) ≈ 0.9096042426`, given as a number because the integral is elliptic.

**Not supported:**
- deg U ≥ 3 (elliptic and beyond; definite integrals give the number);
- on conics without a rational point, residues generating a proper subfield of Q(ρ);
- two radicands in one expression;
- roots of rational functions;
- roots together with exp, ln, sin or cos;
- definite integrals with roots at multiples of π, or over a RootSum (the number only).

**Tests.**
- `test_reg_integrate_sqrt` pins 31 answers (14 by the reduction: powers of U, linear and quadratic poles, a RootSum), each checked by its numeric derivative, plus 7 unsupported inputs and 4 simplify/differentiate forms.
- `test_prop_integrate_sqrt` integrates `r(x)·y + s(x)/y`, `D(p·y + q)`, and on conics without a point `p·y + q/y + k/((x + e)·y)` (sometimes `+ (x + c)/((x² + g)^j·y)`), over four kinds of U, each checked numerically.
- `test_reg_definite_sqrt` pins 25 exact values (9 on infinite ranges), 6 errors (not real, a pole, a non-integrable end, not real on a half-line, two divergent at infinity) and 3 numbers (elliptic on [0, 1] and [0, ∞), a RootSum). `test_prop_definite_sqrt` integrates 60 random integrands on [0, 1/2] (on [−1, 1/2] with a pole (x − h)² + g centred on U for conics without a point, five radicands: c − x², 7 − (x − 1)², 2·(7 − x²)), all exact, and 30 on [0, ∞), at least a third exact (26 of 30 at the default seed; the rest are RootSums from x² + g).
- The SymPy stream `check_definite_sqrt` compares exact values with `mpmath.quad` to 30 digits. On its original 120 cases all are now exact and match; with `1/((x² + g)·y)` terms and shifted or scaled radicands without a point (`6 + 2x − x²`, `14 − 2x²`) added, 85 of 120 are exact; the rest are numeric only, all but 3 from RootSums (those 3 are nested 0/0 atan arguments at a rational point, numeric on `main` too).
- `check_infinite_sqrt` does the same on [c, ∞) and (−∞, ∞), with the reference computed through x = c + 1/v² (and the conjugate form for `(y − √k·x)²`, which loses digits at large x): 115 of 150 exact and matching, 35 numeric (RootSums), 0 failures. Before ln pairs were merged for the limit at ∞ (`−c·ln|A| + c·ln|B|` → `c·ln|A/B|`, each infinite alone), 17 of these were numeric.
- `check_sqrt` differentiates every answer back to 35 digits, or 60 for RootSums at points inside the domain: of 400, 316 elementary, 78 RootSum, and 6 perfect-square radicands correctly refused.

## 12. CLI Pipeline

The `cli/main.c` entry point ties everything together:

```
./polysolver <action> "<expression>"
```

The pipeline for every action:

```
1. parse_and_simplify(input)        (trees with a square root go to algsqrt.c, with exp/ln/sin/cos to risch.c)
   │
   ├─ parse_expression(&input)   ← shunting-yard → Expr* tree
   │
   ├─ simplify(tree)             ← canonicalize → EXPR_POLY or EXPR_RATIONAL_FN
   │
   └─ free_expr(tree)            ← free the intermediate parse tree
   │
   Returns: canonical Expr*

2. Dispatch on action:

   "simplify"      → expr_to_string(expr) → print
   "factor"        → factorize(expr->poly) → factorization_to_string → print
   "differentiate" → differentiate(expr) → expr_to_string → print
   "integrate"     → integrate(expr) → expr_to_string → print
   "solve"         → solve_polynomial(expr->poly) → solutions_to_string → print

3. free_expr(expr) → clean up
```

**Full trace for** `./polysolver factor "x^2 - 5x + 6"`:

```
1. Lexer produces: x, ^, 2, -, 5, x, +, 6
2. Shunting-yard builds: SUB(POW(x, 2), ADD(MUL(5, x), NEG(6)))
   ... actually: SUB(ADD(POW(x,2), MUL(-5, x)), CONST(6))
   ... let's trace it properly via the implicit mul rules:
   
   Actually the canonicalization handles the details. The parser produces
   an Expr tree, and canonicalize() reduces it to:
   
   EXPR_POLY: [1*x^2, -5*x^1, 6*x^0]

3. factorize() is called on this polynomial:
   - Yun: f is square-free → one part, multiplicity 1
   - Berlekamp–Zassenhaus: f ≡ x(x + 1) (mod 3), two linear factors; they lift
     to (x - 3) and (x - 2) mod 3^k, and each divides f exactly over Z
   - sorted by root: {factors: [(x - 2), (x - 3)], multiplicities: [1, 1], content: 1}

4. factorization_to_string: "(x - 2)(x - 3)"

5. Output: (x - 2)(x - 3)
```
