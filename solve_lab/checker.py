#!/usr/bin/env python3
"""Exact integer checker for EQUATIONS.txt.

Transforms each equation `LHS = 0` into a Python expression over an integer
array `v` (v[i] stands for x_i), compiles it once, and evaluates with exact
Python big integers.  No floating point anywhere.

Usage:
  python3 checker.py                       # check all-zeros
  python3 checker.py assignment.json       # check a full/partial assignment
  python3 checker.py --report N            # show first N failing line indices

Assignment JSON: {"x_0": 3, "x_1": -7, ...}.  Missing vars default to 0.
Exit code 0 iff every equation evaluates to exactly 0.
"""
import sys, re, json, time

EQ_PATH = __file__.rsplit('/', 1)[0] + '/../EQUATIONS.txt'
NVARS = 38748
VAR_RE = re.compile(r'x_(\d+)')

def load_equations(path=EQ_PATH):
    """Return (codes, varsets): compiled code objects and per-eq variable id sets."""
    codes = []
    varsets = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            lhs = line.rsplit('=', 1)[0]
            ids = tuple(int(m) for m in VAR_RE.findall(lhs))
            expr = VAR_RE.sub(r'v[\1]', lhs)
            codes.append(compile(expr, '<eq>', 'eval'))
            varsets.append(frozenset(ids))
    return codes, varsets

def evaluate_all(codes, v):
    """Return list of failing indices (equations whose LHS != 0)."""
    ns = {'v': v, '__builtins__': {}}
    fails = []
    for i, c in enumerate(codes):
        if eval(c, ns) != 0:
            fails.append(i)
    return fails

def load_assignment(path):
    with open(path) as f:
        d = json.load(f)
    v = [0] * NVARS
    for k, val in d.items():
        idx = int(k[2:]) if k.startswith('x_') else int(k)
        v[idx] = int(val)
    return v

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    report = 20
    if '--report' in sys.argv:
        report = int(sys.argv[sys.argv.index('--report') + 1])
    t0 = time.time()
    codes, varsets = load_equations()
    print(f"[checker] loaded {len(codes)} equations in {time.time()-t0:.1f}s", file=sys.stderr)
    if args:
        v = load_assignment(args[0])
        src = args[0]
    else:
        v = [0] * NVARS
        src = "all-zeros"
    t1 = time.time()
    fails = evaluate_all(codes, v)
    n = len(codes)
    ok = n - len(fails)
    print(f"[checker] assignment={src}")
    print(f"[checker] satisfied {ok}/{n}  ({len(fails)} failing)  eval={time.time()-t1:.1f}s")
    if fails:
        print(f"[checker] first {min(report,len(fails))} failing line indices: {fails[:report]}")
        print("RESULT: FAIL")
        return 1
    print("RESULT: OK  — all equations satisfied exactly in Z")
    return 0

if __name__ == '__main__':
    sys.exit(main())
