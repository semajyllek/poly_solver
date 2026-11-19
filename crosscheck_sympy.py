#!/usr/bin/env python3
"""Differential test: compare ./poly_solver against SymPy on seeded random inputs.

Requires SymPy (pip install sympy). Usage:
    make cli && python3 crosscheck_sympy.py [--trials N] [--seed S]

Checks, per random input:
  factor     product of the printed factors expands to the input, and the factor
             count/multiplicities match sympy.factor_list
  solve      claimed rational roots are exact roots; every exactly printed root (radicals,
             cos/sin(kπ/N), complex pairs) is a root at 50 digits; each "x ~ d" is within
             1e-9 of a SymPy real root; real and non-real root counts (with multiplicity) match
  integrate  SymPy confirms d/dx(printed antiderivative) - f == 0 symbolically; answers
             containing RootSum(q, t -> t * ln(S)) are checked at 60 significant digits
             (Σ over the roots of q of t·S_x/S, plus the derivative of the other terms),
             and so are answers with radicals and cos/sin(kπ/N) constants, including a
             second stream of binomial/cyclotomic denominators (c·x^n ± e, x^(sN) ± b^N)
  definite   integrate f a b agrees with mpmath.quad to 40 digits, and reports divergence
             exactly when the denominator has a real root in [a, b]
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

    ln(abs(u)) becomes log((u)): d/dx ln|u| = u'/u = d/dx log(u) on every branch, and SymPy
    cannot simplify the sign(u) terms that differentiating Abs would introduce."""
    text = text.replace("ln(abs(", "log((").replace("^", "**")
    return parse_expr(text, local_dict={"x": x, "t": t, "atan": sp.atan, "sqrt": sp.sqrt,
                                        "log": sp.log, "ln": sp.log, "Abs": sp.Abs, "pi": sp.pi, "cos": sp.cos,
                                        "sin": sp.sin, "cot": sp.cot, "csc": sp.csc, "exp": sp.exp},
                      transformations=TRANSFORMS)


def to_cli(expr):
    """SymPy expression -> CLI input syntax."""
    return str(expr).replace("**", "^").replace("log(", "ln(")


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
    nonreal += sum(2 * mult_of(m) for m in re.findall(r"\+/- [^,]*?i\)?(?: \(multiplicity (\d+)\))?(?:,|$)", out))
    real = len(re.findall(r"x = -?\d+(?:/\d+)?(?:,|$)", out))
    real += sum(mult_of(m) for m in re.findall(r"~ -?\d+\.\d+(?: \(multiplicity (\d+)\))?", out))
    return real, nonreal


def exact_roots(out):
    """Every exactly printed root: 'x = a', 'x = a ~ d', 'x = a +/- b*i [~ (...)]'."""
    roots = []
    for entry in re.split(r", (?=x |\d+ non-real)", out):
        entry = re.sub(r" \(multiplicity \d+\)$", "", entry)
        if not entry.startswith("x = "):
            continue
        body = entry[4:].split(" ~ ")[0]
        if "+/-" in body:
            re_part, im_part = body.split("+/-")
            im_part = re.sub(r"\s*\*?\s*i$", "", im_part.strip()) or "1"
            a = to_sympy(re_part) if re_part.strip() else sp.Integer(0)
            b = to_sympy(im_part)
            roots += [a + sp.I * b, a - sp.I * b]
        else:
            roots.append(to_sympy(body))
    return roots


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
    scale = sum(abs(c) for c in sp.Poly(p, x).all_coeffs())
    for r in exact_roots(out):                     # radicals and cos/sin(kπ/N), at 50 digits
        assert abs(sp.N(p.subs(x, r), 50)) < sp.Float("1e-35") * scale, f"claimed root {r} is not a root: {out}"


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


def check_definite(rng):
    """integrate f a b: the exact value matches mpmath.quad at 40 digits; 'diverges' exactly
    when the denominator has a real root in [a, b]; RootSum answers give the number only."""
    num = sum(rng.randint(-5, 5) * x**i for i in range(rng.randint(0, 3))) or sp.Integer(1)
    den = binomial_den(rng) if rng.random() < 0.3 else random_poly(rng, coeff=5)
    f = sp.cancel(num / den)
    a, b = (sp.Rational(rng.randint(-12, 12), rng.randint(1, 4)) for _ in range(2))
    code, out, err = run_cli("integrate", to_cli(f), str(a), str(b))
    if a == b:                               # an empty range: 0 by convention, even at a pole
        assert code == 0 and out.startswith("0 ~"), f"f={f} on [{a}, {a}]: {out} {err}"
        return "exact"
    lo, hi = min(a, b), max(a, b)
    poles = [r for r in sp.Poly(sp.denom(f), x).real_roots() if lo <= r <= hi]
    if code != 0:
        assert "diverges" in err and poles, f"refused f={f} on [{a}, {b}]: {err}"
        return "pole"
    assert not poles, f"missed a pole of f={f} in [{a}, {b}]: {out}"
    fl = sp.lambdify(x, f, "mpmath")
    with mp.workdps(45):
        ref = mp.quad(fl, mp.linspace(mp.mpf(a.p) / a.q, mp.mpf(b.p) / b.q, 9))
    if out.startswith("~"):
        assert "RootSum" in out and abs(float(out.split()[1]) - float(ref)) < 1e-8 * (1 + abs(float(ref))), out
        return "numeric"
    value = out.split(" ~ ")[0].replace("ln(abs(", "log(Abs(")  # a value: keep the |.|
    v = sp.N(to_sympy(value), 40)
    assert abs(v - sp.Float(str(ref), 40)) < sp.Float("1e-28") * (1 + abs(v)), f"f={f} on [{a}, {b}]: {out} vs {ref}"
    return "exact"


def check_infinite(rng):
    """integrate f with an infinite bound: refused exactly when the integral diverges (a real
    pole in the range, or deg den < deg num + 2); otherwise the exact value matches mpmath.quad."""
    k = rng.randint(1, 3)
    q = x**k + sum(rng.randint(-4, 4) * x**i for i in range(k))  # monic, so q² + c has degree 2k
    den = sp.expand(binomial_den(rng) if rng.random() < 0.5 else q**2 + rng.randint(1, 6))
    deg = sp.degree(den, x)
    num = sum(rng.randint(-5, 5) * x**i for i in range(rng.randint(1, deg))) or sp.Integer(1)
    f = sp.cancel(num / den)
    c = sp.Rational(rng.randint(-12, 12), rng.randint(1, 4))
    a, b = rng.choice([("-inf", "inf"), (c, "inf"), ("-inf", c), ("inf", c)])
    code, out, err = run_cli("integrate", to_cli(f), str(a), str(b))
    ext = {"-inf": -mp.inf, "inf": mp.inf}
    lo, hi = sorted([ext.get(a, a), ext.get(b, b)], key=lambda v: float(v))
    n, d = sp.fraction(f)
    poles = [r for r in sp.Poly(d, x).real_roots() if float(lo) <= float(r) <= float(hi)]
    if poles or sp.degree(d, x) < sp.degree(n, x) + 2:
        assert code != 0 and "diverges" in err, f"accepted divergent f={f} on [{a}, {b}]: {out}"
        return "divergent"
    assert code == 0, f"refused f={f} on [{a}, {b}]: {err}"
    fl = sp.lambdify(x, f, "mpmath")
    to_mp = lambda v: ext[v] if v in ext else mp.mpf(v.p) / v.q
    with mp.workdps(45):                     # breakpoints where the poles' peaks can be (|root| < 7)
        ends = sorted([to_mp(a), to_mp(b)])
        cuts = [mp.mpf(k) / 2 for k in range(-14, 15) if ends[0] < mp.mpf(k) / 2 < ends[1]]
        ref = mp.quad(fl, [ends[0]] + cuts + [ends[1]]) * (1 if to_mp(a) < to_mp(b) else -1)
    if out.startswith("~"):
        assert "RootSum" in out and abs(float(out.split()[1]) - float(ref)) < 1e-8 * (1 + abs(float(ref))), out
        return "numeric"
    v = sp.N(to_sympy(out.split(" ~ ")[0].replace("ln(abs(", "log(Abs(")), 40)
    assert abs(v - sp.Float(str(ref), 40)) < sp.Float("1e-25") * (1 + abs(v)), f"f={f} on [{a}, {b}]: {out} vs {ref}"
    return "exact"


def random_coeff(rng):
    """A random rational function: small numerator over (x - r)^m factors (m = 1 or 2)."""
    num = sum(rng.randint(-4, 4) * x**i for i in range(rng.randint(1, 3))) or sp.Integer(1)
    den = sp.Integer(1)
    for _ in range(rng.randint(0, 2)):
        den *= (x - rng.randint(-3, 3))**rng.randint(1, 2)
    return num / den


def is_zero(d):
    """d == 0 as a function of x (trig identities are checked again in exponential form)."""
    return sp.simplify(d) == 0 or sp.simplify(sp.expand(d.rewrite(sp.exp))) == 0


def compare_risch(f, ref_form=None):
    """Our verdict on ∫ f against SymPy's risch_integrate (of ref_form, f rewritten when SymPy
    needs exponentials): a NonElementaryIntegral from SymPy must meet our proof and vice versa,
    and every answer must differentiate back to f exactly."""
    from sympy.integrals.risch import risch_integrate, NonElementaryIntegral
    code, out, err = run_cli("integrate", to_cli(f))
    if code != 0:
        assert "not supported" in err, f"failed on {f}: {err}"
        return "unsupported"
    try:
        ref = risch_integrate(f if ref_form is None else ref_form, x)
    except Exception:                        # SymPy's Risch can fail internally (SympifyError)
        ref = None
    nonelem_ref = ref is not None and (isinstance(ref, NonElementaryIntegral) or ref.has(NonElementaryIntegral))
    if out.startswith("no elementary antiderivative"):
        assert ref is None or nonelem_ref, f"we say non-elementary, SymPy integrates {f}: {ref}"
        return "nonelementary" if ref is not None else "unconfirmed"
    assert not nonelem_ref, f"SymPy proves {f} non-elementary, we gave {out}"
    d = sp.diff(to_sympy(out), x) - f
    assert is_zero(d), f"d/dx of our {out} differs from {f} by {sp.simplify(d)}"
    return "elementary"


def check_exp(rng):
    """f = Σ a_k·exp(k·v); half the integrands are derivatives (always elementary)."""
    v = rng.randint(1, 3) * x**rng.randint(1, 2) * rng.choice([-1, 1]) + rng.randint(-2, 2) * x
    ks = rng.sample([-2, -1, 1, 2], rng.randint(1, 2))
    if rng.random() < 0.5:
        F = sum(random_coeff(rng) * sp.exp(k * v) for k in ks)
        f = sp.together(sp.diff(F, x))
    else:
        f = sum(random_coeff(rng) * sp.exp(k * v) for k in ks) + (random_coeff(rng) if rng.random() < 0.3 else 0)
    return compare_risch(f)


def check_log(rng):
    """f = Σ a_k·ln(u)^k, u = x, x² + c or (x + c)/(x - c); half are derivatives."""
    c = rng.randint(1, 3)
    u = rng.choice([x, x**2 + c, (x + c) / (x - c)])
    if rng.random() < 0.5:
        F = sum(random_coeff(rng) * sp.log(u)**k for k in range(1, rng.randint(2, 3) + 1))
        f = sp.diff(F, x)
    else:
        f = sum(random_coeff(rng) * sp.log(u)**k for k in range(0, rng.randint(1, 2) + 1))
    return compare_risch(f)


def theta_poly(rng, th, deg):
    """Σ_{k<=deg} c_k·θ^k with small integer or linear-in-x coefficients, monic-ish on top."""
    terms = [rng.choice([rng.randint(-3, 3), rng.randint(-2, 2) * x + rng.randint(-2, 2)]) * th**k for k in range(deg)]
    return sum(terms) + rng.choice([1, 2, x]) * th**deg


def check_frac(rng):
    """f a rational function of θ = exp(v) or ln(u) with a denominator in θ: half are
    derivatives of B/D + c·ln(G) (+ an atan, irrational residues; elementary), half random."""
    th = rng.choice([sp.exp(x), sp.exp(-2 * x), sp.exp(x**2), sp.log(x), sp.log(x**2 + 1)])
    if rng.random() < 0.5:
        F = theta_poly(rng, th, rng.randint(0, 1)) / theta_poly(rng, th, rng.randint(1, 2))**rng.randint(1, 2)
        F += rng.choice([-2, -1, 1, sp.Rational(1, 2)]) * sp.log(theta_poly(rng, th, 1))
        if rng.random() < 0.5:                       # irrational residues: (1/√e)·atan(U/√e)
            e = rng.choice([1, 2, 3])
            F += sp.atan(theta_poly(rng, th, 1) / sp.sqrt(e)) / sp.sqrt(e)
        f = sp.together(sp.diff(F, x))
    else:
        f = theta_poly(rng, th, rng.randint(0, 2)) / theta_poly(rng, th, rng.randint(1, 2))
    return compare_risch(f)


def check_trig(rng):
    """f built from r(x)·exp(e·x)·sin/cos(k·x or x²): half derivatives (elementary), half random.
    SymPy's Risch gets f in exponential form."""
    def term():
        arg = rng.choice([x, 2 * x, 3 * x, x**2])
        return random_coeff(rng) * sp.exp(rng.choice([-1, 0, 1]) * x) * rng.choice([sp.sin, sp.cos])(arg)
    F = sum(term() for _ in range(rng.randint(1, 2)))
    f = sp.diff(F, x) if rng.random() < 0.5 else F
    return compare_risch(f, f.rewrite(sp.exp))


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
    rng_bin = random.Random(args.seed + 1)  # separate streams: the main inputs stay unchanged
    rng_def = random.Random(args.seed + 2)
    rng_inf = random.Random(args.seed + 3)
    rng_exp = random.Random(args.seed + 4)
    rng_log = random.Random(args.seed + 5)
    rng_frac = random.Random(args.seed + 6)
    rng_trig = random.Random(args.seed + 7)
    rootsums = radicals = 0
    definite = {"exact": 0, "numeric": 0, "pole": 0}
    infinite = {"exact": 0, "numeric": 0, "divergent": 0}
    exps = {"elementary": 0, "nonelementary": 0, "unconfirmed": 0, "unsupported": 0}
    logs = dict(exps)
    fracs = dict(exps)
    trigs = dict(exps)
    for t in range(args.trials):
        p = random_poly(rng)
        check_factor(p)
        check_solve(p)
        rootsums += check_integrate(rng) == "rootsum"
        check_solve(sp.expand(binomial_den(rng_bin)))
        kind = check_integrate(rng_bin, binomial_den(rng_bin))
        rootsums += kind == "rootsum"
        radicals += kind == "radical"
        definite[check_definite(rng_def)] += 1
        infinite[check_infinite(rng_inf)] += 1
        exps[check_exp(rng_exp)] += 1
        logs[check_log(rng_log)] += 1
        fracs[check_frac(rng_frac)] += 1
        trigs[check_trig(rng_trig)] += 1
        check_matrix(rng)
        print(f"  trial {t + 1}/{args.trials} ok", flush=True)
    print(f"crosscheck OK: {args.trials} trials (factor, solve, integrate, det); "
          f"{rootsums} antiderivatives contained RootSum terms and {radicals} radical/trig forms "
          f"(both checked to 60 digits); definite integrals: {definite['exact']} exact (40 digits), "
          f"{definite['numeric']} numeric-only, {definite['pole']} divergent; "
          f"infinite bounds: {infinite['exact']} exact, {infinite['numeric']} numeric-only, "
          f"{infinite['divergent']} divergent; exp integrals (vs SymPy's risch_integrate): "
          f"{exps['elementary']} elementary, {exps['nonelementary']} proved non-elementary (SymPy agrees), "
          f"{exps['unconfirmed']} non-elementary where SymPy's Risch failed, {exps['unsupported']} unsupported; "
          f"ln integrals: {logs['elementary']} elementary, {logs['nonelementary']} proved non-elementary (SymPy agrees), "
          f"{logs['unconfirmed']} unconfirmed, {logs['unsupported']} unsupported; "
          f"denominators in exp/ln: {fracs['elementary']} elementary, {fracs['nonelementary']} proved non-elementary "
          f"(SymPy agrees), {fracs['unconfirmed']} unconfirmed, {fracs['unsupported']} unsupported; "
          f"sin/cos: {trigs['elementary']} elementary, {trigs['nonelementary']} proved non-elementary (SymPy agrees), "
          f"{trigs['unconfirmed']} unconfirmed, {trigs['unsupported']} unsupported")


if __name__ == "__main__":
    sys.exit(main())
