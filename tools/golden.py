#!/usr/bin/env python3
"""Golden corpus: exact CLI outputs, replayed in seconds as the fast regression lane.

    make cli && python3 golden.py                 # replay golden.jsonl: same exit code and output
    python3 golden.py --record --seeds 1 2 3      # rebuild it from crosscheck_sympy's inputs (needs sympy)

Replay needs no SymPy: each line is one CLI call with the exit code and stripped stdout it
gave when recorded. Re-record only when an output change is intended (and checked).
"""
import argparse
import json
import resource
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLI = [os.path.join(ROOT, "polysolver")]  # the binary under test (--cli)
CORPUS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden.jsonl")


def run(args):
    try:
        r = subprocess.run([CLI[0], *args], capture_output=True, text=True, timeout=60)
        return r.returncode, r.stdout.strip()
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT"


def replay(path, jobs):
    calls = [json.loads(line) for line in open(path)]
    bad = 0
    with ThreadPoolExecutor(jobs) as ex:
        for c, (code, out) in zip(calls, ex.map(lambda c: run(c["args"]), calls)):
            if (code, out) != (c["code"], c["out"]):
                bad += 1
                if bad <= 5:
                    print("DIFF", c["args"], "\n  was:", c["out"][:400], "\n  now:", out[:400])
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    print(f"golden: {len(calls)} calls, {bad} differences, CPU {ru.ru_utime + ru.ru_stime:.1f} s")
    return bad == 0


def record(path, seeds, trials):
    import crosscheck_sympy as cs  # the crosscheck's generators pick the inputs
    seen, rows = set(), []
    orig = cs.run_cli

    def rec(*args):
        r = orig(*args)
        if args not in seen:  # one row per distinct call
            seen.add(args)
            rows.append({"args": list(args), "code": r[0], "out": r[1]})
        return r

    cs.run_cli = rec
    for s in seeds:
        sys.argv = ["crosscheck_sympy.py", "--trials", str(trials), "--seed", str(s)]
        cs.main()
    with open(path, "w") as f:
        for row in rows:
            f.write(json.dumps(row) + "\n")
    print(f"golden: recorded {len(rows)} calls")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", action="store_true")
    ap.add_argument("--cli", help="the polysolver binary to replay against (default: ./polysolver)")
    ap.add_argument("--seeds", type=int, nargs="*", default=[1, 2, 3])
    ap.add_argument("--trials", type=int, default=40)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    if a.cli:
        CLI[0] = os.path.abspath(a.cli)
    if a.record:
        record(CORPUS, a.seeds, a.trials)
    else:
        sys.exit(0 if replay(CORPUS, a.jobs) else 1)


if __name__ == "__main__":
    main()
