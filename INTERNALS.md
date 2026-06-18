# poly_solver Internals

A comprehensive guide to how every algorithm, data structure, and module in poly_solver works, illustrated with toy examples.

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

poly_solver is a symbolic algebra engine written in C. It takes a string like `"x^2 - 5x + 6"`, parses it into an internal tree representation, reduces it to a canonical polynomial form, and then applies operations—simplification, factoring, differentiation, integration, or equation solving.

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
rational_num.c/h          ← standalone, no dependencies
    │
    ▼
polynomial.c/h            ← depends on rational_num
    │
    ▼
expression.c/h            ← depends on polynomial, rational_num
    │
    ├──► factor.c/h       ← depends on polynomial, expression
    │
    ▼
parser.c/h                ← depends on expression, token
    │
    ▼
solver.c/h                ← depends on polynomial, factor
    │
    ▼
main.c                    ← depends on everything, wires CLI
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

The `Rational` type is a **value type** (16 bytes, two `long long` fields) passed by value on the stack. It never needs freeing. The only heap-allocated rational output is `rat_to_string()`, which the caller must free.

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
typedef long long rat_int_t;  // typedef allows future GMP swap

typedef struct {
    rat_int_t num;   // numerator (carries the sign)
    rat_int_t den;   // denominator (always > 0 after normalization)
} Rational;
```

Every `Rational` is kept in **lowest terms** with a **positive denominator**. This normalization is enforced by the `normalize()` function called by every constructor and arithmetic operation.

### 4.3 Normalization via Euclidean GCD

The Euclidean algorithm computes `gcd(a, b)` by repeatedly replacing the larger number with the remainder of dividing the two:

```
gcd(12, 8):
    12 mod 8 = 4   → gcd(8, 4)
     8 mod 4 = 0   → gcd(4, 0) → answer: 4
```

In code (`gcd_int`):
```c
rat_int_t gcd_int(rat_int_t a, rat_int_t b) {
    if (a < 0) a = -a;  // work with absolute values
    if (b < 0) b = -b;
    while (b != 0) {
        rat_int_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}
```

Normalization then divides both numerator and denominator by their GCD and ensures the denominator is positive:

```
normalize(6, -4):
    den < 0 → flip signs → (-6, 4)
    gcd(6, 4) = 2
    result: -6/2 = -3,  4/2 = 2  → Rational{-3, 2}
```

### 4.4 Arithmetic with overflow mitigation

**Addition** uses the identity `a/b + c/d = (a*(d/g) + c*(b/g)) / (b/g*d)` where `g = gcd(b, d)`. By dividing out the GCD of the denominators first, intermediate products stay smaller:

```
1/6 + 1/4:
    g = gcd(6, 4) = 2
    lcm_den = 6/2 * 4 = 12
    num = 1*(4/2) + 1*(6/2) = 2 + 3 = 5
    normalize(5, 12) → 5/12
```

**Multiplication** cross-reduces before multiplying: `(a/b) * (c/d)` first computes `g1 = gcd(a, d)` and `g2 = gcd(c, b)`, then multiplies the reduced numerators and denominators:

```
(2/3) * (3/4):
    g1 = gcd(2, 4) = 2     g2 = gcd(3, 3) = 3
    num = (2/2) * (3/3) = 1 * 1 = 1
    den = (3/3) * (4/2) = 1 * 2 = 2
    result: 1/2
```

Without cross-reduction, the intermediate would be `6/12` — same result after normalization, but with larger intermediates that risk overflow.

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

Addition creates a new polynomial by inserting all terms from both operands. Because `add_term` handles sorting and like-term combination, the algorithm is simple:

```c
Polynomial* poly_add(poly1, poly2) {
    result = new empty polynomial
    for each term in poly1: add_term(result, term.coeff, term.exp)
    for each term in poly2: add_term(result, term.coeff, term.exp)
    return result
}
```

**Example:** `(3x^2 + 4x) + (5x^2 + 1)`:

```
Step 1: Insert terms from poly1:
    add 3x^2 → [3, exp=2]
    add 4x   → [3, exp=2] → [4, exp=1]

Step 2: Insert terms from poly2:
    add 5x^2 → found exp=2, combine: 3+5=8 → [8, exp=2] → [4, exp=1]
    add 1    → insert at end → [8, exp=2] → [4, exp=1] → [1, exp=0]

Result: 8x^2 + 4x + 1
```

Subtraction works identically but negates the coefficients from the second polynomial.

### 5.4 Polynomial multiplication

Uses the distributive property: multiply every term of poly1 by every term of poly2, then collect via `add_term`:

```
(2x + 3) * (x + 1):

    2x * x = 2x^2       → add_term(result, 2, 2)
    2x * 1 = 2x         → add_term(result, 2, 1)
    3  * x = 3x          → add_term(result, 3, 1)  ← combines with existing: 2+3=5
    3  * 1 = 3           → add_term(result, 3, 0)

Result: 2x^2 + 5x + 3
```

Time complexity: O(n*m) where n, m are the number of terms.

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

The convenience wrappers `poly_divide` and `poly_mod` call `poly_divmod` and discard the unwanted part.

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
Level 4: unary minus (-)          right-associative
Level 3: power (^)                right-associative
Level 2: multiply (*), divide (/) left-associative
Level 1: add (+), subtract (-)    left-associative
Level 0: parentheses, sentinel
```

Power is right-associative so `x^2^3` parses as `x^(2^3)`, not `(x^2)^3`.

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

### 9.1 Overview of the factoring pipeline

```
Input polynomial
    │
    ▼
┌────────────────────┐
│ Extract content    │  GCD of all coefficients
│ (poly_content)     │  → separate leading scalar
└────────────────────┘
    │  primitive polynomial (content = 1)
    ▼
┌────────────────────┐
│ Clear denominators │  multiply so all coefficients are integers
└────────────────────┘
    │  integer-coefficient polynomial
    ▼
┌────────────────────┐
│ Rational Root      │  test p/q candidates where p | constant, q | leading
│ Theorem            │
└────────────────────┘
    │  list of rational roots
    ▼
┌────────────────────┐
│ Synthetic Division │  divide out each (x - root) factor
└────────────────────┘
    │  factored form + irreducible remainder
    ▼
  Factorization result
```

### 9.2 Content extraction

The **content** of a polynomial is the GCD of all its coefficients. For rational coefficients, the GCD is computed as:

```
gcd(a/b, c/d) = gcd(a, c) / lcm(b, d)
```

**Example:** Content of `(1/2)x^2 - (3/4)x + (1/4)`:

```
Coefficients: 1/2, -3/4, 1/4
Absolute values: 1/2, 3/4, 1/4

gcd(1/2, 3/4):
    gcd(1, 3) = 1
    lcm(2, 4) = 4
    → 1/4

gcd(1/4, 1/4):
    gcd(1, 1) = 1
    lcm(4, 4) = 4
    → 1/4

Content = 1/4
```

Dividing all coefficients by `1/4` gives the primitive part: `2x^2 - 3x + 1`.

### 9.3 The Rational Root Theorem

For a polynomial `a_n*x^n + ... + a_1*x + a_0` with integer coefficients, any rational root `p/q` (in lowest terms) must satisfy:

- `p` divides the constant term `a_0`
- `q` divides the leading coefficient `a_n`

**Example:** Finding roots of `2x^3 - 3x^2 - 8x - 3`:

```
a_0 = -3, a_n = 2

Divisors of |-3|: 1, 3
Divisors of |2|:  1, 2

Candidate roots (±p/q):
    ±1/1, ±3/1, ±1/2, ±3/2

Test each by evaluating the polynomial:
    f(1)   = 2 - 3 - 8 - 3 = -12  ✗
    f(-1)  = -2 - 3 + 8 - 3 = 0   ✓  root found!
    
After dividing out (x + 1), test remaining polynomial 2x^2 - 5x - 3:
    f(3)   = 18 - 15 - 3 = 0      ✓  root found!
    
After dividing out (x - 3), remainder: 2x + 1
    f(-1/2) = -1 + 1 = 0          ✓  root found!
```

The implementation clears denominators first (multiplying all coefficients by the LCM of their denominators) so that the rational root theorem applies to integer coefficients.

When a constant term is zero, the algorithm immediately records `x = 0` as a root and divides out `x` before proceeding.

When a root is found, the algorithm re-tests the same candidate to detect **repeated roots**: if `x = r` is a root of `f(x)`, it divides out `(x - r)` and checks if `r` is still a root of the quotient.

### 9.4 Synthetic division

Synthetic division is an efficient algorithm for dividing a polynomial by `(x - r)`. It's much faster than full polynomial long division for this special case.

**Example:** Dividing `2x^3 - 3x^2 - 8x - 3` by `(x - 3)`:

```
Write coefficients in descending order: [2, -3, -8, -3]
Root: r = 3

    Step 0: bring down first coefficient:    result[0] = 2
    Step 1: result[1] = -3 + 2*3 = -3 + 6 = 3
    Step 2: result[2] = -8 + 3*3 = -8 + 9 = 1

Quotient coefficients: [2, 3, 1] → 2x^2 + 3x + 1
(The remainder would be -3 + 1*3 = 0, confirming 3 is a root)
```

The implementation converts the sparse linked list to a dense coefficient array for the synthetic division, then converts back.

### 9.5 The `Factorization` result

```c
typedef struct {
    Polynomial** factors;     // array of irreducible factors
    int* multiplicities;      // multiplicity of each factor
    int count;                // number of distinct factors
    Rational content;         // leading scalar (GCD of original coefficients)
} Factorization;
```

**Example:** Factoring `6x^2 + 12x + 6`:

```
Content: 6
Primitive part: x^2 + 2x + 1

Rational roots of x^2 + 2x + 1:
    f(-1) = 1 - 2 + 1 = 0 → root -1
    Divide by (x + 1): quotient = x + 1
    f(-1) again: -1 + 1 = 0 → repeated root!
    Divide by (x + 1): quotient = 1 (constant)

Result:
    content = 6
    factors = [(x + 1)]
    multiplicities = [2]

Displayed as: 6(x + 1)^2
```

### 9.6 Handling irreducible remainders

After extracting all rational roots, the remaining polynomial may be irreducible over Q (has no rational roots). It is added as a factor with multiplicity 1.

**Example:** Factoring `x^6 - 1`:

```
Rational roots: 1, -1

After dividing out (x - 1) and (x + 1):
    x^6 - 1 / (x - 1) = x^5 + x^4 + x^3 + x^2 + x + 1
    that / (x + 1) = x^4 + x^2 + 1

x^4 + x^2 + 1 has no rational roots → kept as irreducible factor.

Result: (x - 1)(x + 1)(x^4 + x^2 + 1)
```

(Mathematically, `x^4 + x^2 + 1 = (x^2 + x + 1)(x^2 - x + 1)`, but neither quadratic has rational roots, so the factorizer cannot split it further without more advanced algorithms.)

---

## 10. Equation Solving

### 10.1 Strategy: factor first, then solve each factor

The solver delegates to the factoring engine, then solves each irreducible factor individually based on its degree:

```
solve_polynomial(p):
    f = factorize(p)
    for each factor in f:
        if degree == 1: solve linear
        if degree == 2: solve quadratic (formula)
        if degree > 2:  report as irreducible (no closed-form solution)
```

### 10.2 Linear solving

For `ax + b = 0`:  `x = -b/a` (exact rational result).

### 10.3 Quadratic solving

For `ax^2 + bx + c = 0`, compute the discriminant `D = b^2 - 4ac` using exact rational arithmetic.

Three cases:

**D = 0:** One repeated root: `x = -b/(2a)`

```
x^2 - 2x + 1 = 0:
    D = 4 - 4 = 0
    x = 2/2 = 1 (multiplicity 2)
```

**D > 0, perfect square:** Two exact rational roots via the quadratic formula.

```
x^2 - 5x + 6 = 0:
    D = 25 - 24 = 1 → √1 = 1
    x = (5 ± 1) / 2 → x = 3 or x = 2
```

The `is_perfect_square_rat` function checks if both the numerator and denominator of `D` are perfect squares by computing their floating-point square roots, rounding, and verifying `round(√n)^2 == n`.

**D > 0, not perfect square:** Two irrational roots, returned as doubles.

```
x^2 - 2 = 0:
    D = 0 + 8 = 8 → not a perfect square
    x ≈ (0 ± 2.828) / 2 → x ≈ ±1.414
```

**D < 0:** No real roots. The irreducible quadratic is recorded for display.

```
x^2 + 1 = 0:
    D = 0 - 4 = -4 < 0
    → "(x^2 + 1) = 0 (no rational roots)"
```

### 10.4 Worked example: solving a cubic

Input: `x^3 - 6x^2 + 11x - 6 = 0`

```
Step 1: Factorize
    Content: 1 (monic)
    Rational root theorem candidates: ±1, ±2, ±3, ±6
    
    f(1) = 1 - 6 + 11 - 6 = 0 → root 1
    Synthetic divide by (x - 1): x^2 - 5x + 6
    
    f(2) = 4 - 10 + 6 = 0 → root 2
    Synthetic divide by (x - 2): x - 3
    
    f(3) = 3 - 3 = 0 → root 3
    
    Factorization: (x - 1)(x - 2)(x - 3)

Step 2: Solve each linear factor
    x - 1 = 0 → x = 1
    x - 2 = 0 → x = 2
    x - 3 = 0 → x = 3

Output: "x = 1, x = 2, x = 3"
```

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

### 11.2 Integration

`integrate()` canonicalizes, then delegates to `poly_integral` for polynomial inputs. Rational function integration (partial fractions) is not yet implemented.

---

## 12. CLI Pipeline

The `main.c` entry point ties everything together:

```
./poly_solver <action> "<expression>"
```

The pipeline for every action:

```
1. parse_and_simplify(input)
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

**Full trace for** `./poly_solver factor "x^2 - 5x + 6"`:

```
1. Lexer produces: x, ^, 2, -, 5, x, +, 6
2. Shunting-yard builds: SUB(POW(x, 2), ADD(MUL(5, x), NEG(6)))
   ... actually: SUB(ADD(POW(x,2), MUL(-5, x)), CONST(6))
   ... let's trace it properly via the implicit mul rules:
   
   Actually the canonicalization handles the details. The parser produces
   an Expr tree, and canonicalize() reduces it to:
   
   EXPR_POLY: [1*x^2, -5*x^1, 6*x^0]

3. factorize() is called on this polynomial:
   - Content: 1 (it's monic)
   - Primitive part: x^2 - 5x + 6
   - Rational root candidates: ±1, ±2, ±3, ±6
   - f(1) = 1-5+6 = 2 ≠ 0
   - f(-1) = 1+5+6 = 12 ≠ 0
   - f(2) = 4-10+6 = 0 ✓ → root found
   - Synthetic divide: (x^2-5x+6)/(x-2) = x-3
   - f(3) on (x-3): 3-3=0 ✓ → root found
   - Factorization: {factors: [(x-2), (x-3)], multiplicities: [1, 1], content: 1}

4. factorization_to_string: "(x - 2)(x - 3)"

5. Output: (x - 2)(x - 3)
```
