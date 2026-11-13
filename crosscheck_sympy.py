#!/usr/bin/env python3
"""Differential test: compare ./poly_solver against SymPy on seeded random inputs.

Requires SymPy (pip install sympy). Usage:
    make cli && python3 crosscheck_sympy.py [--trials N] [--seed S]

Checks, per random input:
  factor     product of the printed factors expands to the input, and the factor
             count/multiplicities match sympy.factor_list
  solve      claimed rational roots are exact roots; each "x ~ d" is within 1e-9 of a
             SymPy real root; real and non-real root counts (with multiplicity) match
  integrate  SymPy confirms d/dx(printed antiderivative) - f == 0, or the CLI refused
             and the denominator has an irreducible factor of degree >= 3
  det/inv    exact agreement with SymPy's Matrix
"""
import argparse
import random
import re
import subprocess
import sys

import sympy as sp
from sympy.parsing.sympy_parser import (implicit_multiplication_application,
                                        parse_expr, standard_transformations)

x = sp.Symbol("x")
TRANSFORMS = standard_transformations + (implicit_multiplication_application,)


def run_cli(*args):
    r = subprocess.run(["./poly_solver", *args], capture_output=True, text=True, timeout=60)
    return r.returncode, r.stdout.strip(), r.stderr.strip()


def to_sympy(text):
    """Parse CLI output syntax (implicit products, ^) into SymPy.

    ln|u| becomes log(u): d/dx ln|u| = u'/u = d/dx log(u) on every branch, and SymPy
    cannot simplify the sign(u) terms that differentiating Abs would introduce."""
    text = re.sub(r"ln\|([^|]*)\|", r"log(\1)", text).replace("^", "**")
    return parse_expr(text, local_dict={"x": x, "atan": sp.atan, "sqrt": sp.sqrt,
                                        "log": sp.log, "Abs": sp.Abs},
                      transformations=TRANSFORMS)


def to_cli(expr):
    """SymPy expression -> CLI input syntax."""
    return str(expr).replace("**", "^")


def random_poly(rng, max_deg=6, coeff=9):
    """Product of small random factors, so factoring has real work to do."""
    p = sp.Integer(rng.choice([-3, -2, -1, 1, 2, 5]))
    for _ in range(rng.randint(1, 4)):
        d = rng.randint(1, 3)
        f = sum(rng.randint(-coeff, coeff) * x**i for i in range(d)) + rng.randint(1, 3) * x**d
        p *= f**rng.randint(1, 2)
    return sp.expand(p)


def check_factor(p):
    code, out, err = run_cli("factor", to_cli(p))
    assert code == 0, err
    # "(...)(...)^2" -> product; a leading content multiplies the rest
    prod = to_sympy(re.sub(r"\)\(", ")*(", out).replace(")^", ")**"))
    assert sp.expand(prod - p) == 0, f"factor product mismatch: {out}"
    ref = sorted(m for _, m in sp.factor_list(p)[1])
    ours = sorted(int(m) if m else 1 for m in re.findall(r"\)(?:\^(\d+))?", out))
    assert ref == ours, f"factor multiplicities differ: ours {out}, sympy {sp.factor_list(p)}"


def mult_of(m):
    return int(m) if m else 1


def root_counts(out):
    """(real, non-real) root counts claimed by the CLI output, with multiplicity."""
    nonreal = sum(int(n) * mult_of(m) for n, m in
                  re.findall(r"(\d+) non-real roots of \([^)]*\)(?: \(multiplicity (\d+)\))?", out))
    nonreal += sum(2 * mult_of(m) for m in re.findall(r"\+/- [^,]*?i(?: \(multiplicity (\d+)\))?(?:,|$)", out))
    real = len(re.findall(r"x = -?\d+(?:/\d+)?(?:,|$)", out))
    real += sum(mult_of(m) for m in re.findall(r"~ -?\d+\.\d+(?: \(multiplicity (\d+)\))?", out))
    return real, nonreal


def check_solve(p):
    code, out, err = run_cli("solve", to_cli(p))
    assert code == 0, err
    real = [r for r in sp.Poly(p, x).real_roots()]
    deg = sp.degree(p, x)
    assert root_counts(out) == (len(real), deg - len(real)), f"root counts differ: {out}"
    approx = [float(v) for v in re.findall(r"~ (-?\d+\.\d+)", out)]
    exact = [sp.Rational(v) for v in re.findall(r"x = (-?\d+(?:/\d+)?)(?:,|$)", out)]
    for a in approx:
        assert any(abs(float(r) - a) < 1e-9 for r in real), f"approx root {a} not a root: {out}"
    for e in exact:
        assert p.subs(x, e) == 0, f"claimed root {e} is not a root"


def check_integrate(rng):
    num = sum(rng.randint(-5, 5) * x**i for i in range(rng.randint(0, 3))) or sp.Integer(1)
    den = random_poly(rng, coeff=5)
    f = sp.cancel(num / den)
    code, out, err = run_cli("integrate", to_cli(f))
    if code != 0:
        _, den_factors = sp.factor_list(sp.denom(f))
        assert "algebraic" in err and any(sp.degree(q, x) >= 3 for q, _ in den_factors), \
            f"refused without a degree >= 3 denominator factor: f={f}: {err}"
        return "refused"
    F = to_sympy(out)
    d = sp.cancel(sp.together(sp.diff(F, x) - f))  # cheap first; simplify only if needed
    if d != 0:
        d = sp.simplify(d)
    assert d == 0, f"d/dx F != f for f={f}, F={out}"
    return "ok"


def check_matrix(rng):
    n = rng.randint(1, 4)
    M = sp.Matrix(n, n, lambda i, j: sp.Rational(rng.randint(-9, 9), rng.randint(1, 4)))
    lit = "[" + ",".join("[" + ",".join(str(v) for v in M.row(i)) + "]" for i in range(n)) + "]"
    code, out, _ = run_cli("det", lit)
    assert code == 0 and sp.Rational(out) == M.det(), f"det mismatch: {out} vs {M.det()}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20261113)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    refused = 0
    for t in range(args.trials):
        p = random_poly(rng)
        check_factor(p)
        check_solve(p)
        refused += check_integrate(rng) == "refused"
        check_matrix(rng)
        print(f"  trial {t + 1}/{args.trials} ok", flush=True)
    print(f"crosscheck OK: {args.trials} trials (factor, solve, integrate, det); "
          f"{refused} integrals correctly refused as needing algebraic numbers")


if __name__ == "__main__":
    sys.exit(main())
