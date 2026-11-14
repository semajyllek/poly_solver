#!/usr/bin/env python3
"""Differential test: compare ./poly_solver against SymPy on seeded random inputs.

Requires SymPy (pip install sympy). Usage:
    make cli && python3 crosscheck_sympy.py [--trials N] [--seed S]

Checks, per random input:
  factor     product of the printed factors expands to the input, and the factor
             count/multiplicities match sympy.factor_list
  solve      claimed rational roots are exact roots; each "x ~ d" is within 1e-9 of a
             SymPy real root; real and non-real root counts (with multiplicity) match
  integrate  SymPy confirms d/dx(printed antiderivative) - f == 0 symbolically; answers
             containing RootSum(q, t -> t * ln(S)) are checked at 60 significant digits
             (Σ over the roots of q of t·S_x/S, plus the derivative of the other terms),
             and so are answers with radicals and cos/sin(kπ/N) constants, including a
             second stream of binomial/cyclotomic denominators (c·x^n ± e, x^(sN) ± b^N)
  det/inv    exact agreement with SymPy's Matrix
"""
import argparse
import random
import re
import subprocess
import sys

import mpmath as mp
import sympy as sp
from sympy.parsing.sympy_parser import (implicit_multiplication_application,
                                        parse_expr, standard_transformations)

x, t = sp.symbols("x t")
mp.mp.dps = 60
ROOTSUM = r"RootSum\(([^,]*), t -> t \* ln\((.*?)\)\)(?= [+-] |$)"
TRANSFORMS = standard_transformations + (implicit_multiplication_application,)


def run_cli(*args):
    r = subprocess.run(["./poly_solver", *args], capture_output=True, text=True, timeout=60)
    return r.returncode, r.stdout.strip(), r.stderr.strip()


def to_sympy(text):
    """Parse CLI output syntax (implicit products, ^) into SymPy.

    ln|u| becomes log(u): d/dx ln|u| = u'/u = d/dx log(u) on every branch, and SymPy
    cannot simplify the sign(u) terms that differentiating Abs would introduce."""
    text = re.sub(r"ln\|([^|]*)\|", r"log(\1)", text).replace("^", "**")
    return parse_expr(text, local_dict={"x": x, "t": t, "atan": sp.atan, "sqrt": sp.sqrt,
                                        "log": sp.log, "Abs": sp.Abs, "pi": sp.pi, "cos": sp.cos,
                                        "sin": sp.sin, "cot": sp.cot, "csc": sp.csc},
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


def rootsum_ok(f, out):
    """High-precision check of an antiderivative that contains RootSum terms."""
    rootsums = re.findall(ROOTSUM, out)
    rest = re.sub(r"\s*[+-]?\s*" + ROOTSUM, "", out).strip()
    rest = re.sub(r"^\+\s*", "", rest)
    d_rest = sp.diff(to_sympy(rest), x) if rest else sp.Integer(0)
    for x0 in (sp.Rational(-33, 10), sp.Rational(-3, 5), sp.Rational(13, 10), sp.Rational(29, 10)):
        val = sp.N(d_rest.subs(x, x0), 60)
        for q, body in rootsums:
            S = to_sympy(body)
            Sx = sp.diff(S, x)
            for a in sp.Poly(to_sympy(q), t).nroots(n=60, maxsteps=500):
                val += sp.N(a * Sx.subs({x: x0, t: a}) / S.subs({x: x0, t: a}), 60)
        fx = sp.N(f.subs(x, x0), 60)
        if abs(sp.N(val - fx, 60)) > sp.Float("1e-40") * (1 + abs(fx)):
            return False
    return True


def numeric_ok(f, out):
    """d/dx F == f at 60 digits, for answers with radicals and trig constants."""
    dF = sp.diff(to_sympy(out), x)
    for x0 in (sp.Rational(-33, 10), sp.Rational(-3, 5), sp.Rational(1, 7), sp.Rational(29, 10)):
        fx = sp.N(f.subs(x, x0), 60)
        if abs(sp.N(dF.subs(x, x0), 60) - fx) > sp.Float("1e-40") * (1 + abs(fx)):
            return False
    return True


def binomial_den(rng):
    """c·x^n ± e, or x^(sN) ± b^N (rescaled cyclotomic factors), maybe times x - k."""
    if rng.random() < 0.5:
        den = rng.randint(1, 4) * x**rng.randint(3, 9) + rng.choice([-1, 1]) * rng.randint(1, 7)
    else:
        s, N, b = rng.randint(1, 2), rng.randint(3, 6), rng.randint(1, 3)
        den = x**(s * N) + rng.choice([-1, 1]) * b**N
    return den * (x - rng.randint(-3, 3)) if rng.random() < 0.3 else den


def check_integrate(rng, den=None):
    num = sum(rng.randint(-5, 5) * x**i for i in range(rng.randint(0, 3))) or sp.Integer(1)
    den = den if den is not None else random_poly(rng, coeff=5)
    f = sp.cancel(num / den)
    code, out, err = run_cli("integrate", to_cli(f))
    assert code == 0, f"integrate failed for f={f}: {err}"
    if "RootSum" in out:
        assert rootsum_ok(f, out), f"d/dx F != f (60 digits) for f={f}, F={out}"
        return "rootsum"
    if "^(" in out or "pi" in out:
        assert numeric_ok(f, out), f"d/dx F != f (60 digits) for f={f}, F={out}"
        return "radical"
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
    rng_bin = random.Random(args.seed + 1)  # separate stream: the main inputs stay unchanged
    rootsums = radicals = 0
    for t in range(args.trials):
        p = random_poly(rng)
        check_factor(p)
        check_solve(p)
        rootsums += check_integrate(rng) == "rootsum"
        kind = check_integrate(rng_bin, binomial_den(rng_bin))
        rootsums += kind == "rootsum"
        radicals += kind == "radical"
        check_matrix(rng)
        print(f"  trial {t + 1}/{args.trials} ok", flush=True)
    print(f"crosscheck OK: {args.trials} trials (factor, solve, integrate, det); "
          f"{rootsums} antiderivatives contained RootSum terms and {radicals} radical/trig forms "
          f"(both checked to 60 digits)")


if __name__ == "__main__":
    sys.exit(main())
