"""Simultaneous mod-p solve: drive the bad checks to 0 while HOLDING the solved residuals at 0.
   Controls: every free input.  Iterate: solve -> close -> re-solve."""
import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep, tri7
from gfp import gauss_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
FREE = [u for u in range(L.NVARS) if L.definer.get(u) is None]
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
BITSET = set(int(b) for b in LD)
NONBIT = [u for u in FREE if u not in BITSET]

HOLD = [lambda v: v[3719], lambda v: v[25118],
        lambda v: v[16441] - v[4920], lambda v: v[28955] - v[11408],
        lambda v: v[2751] - v[1085], lambda v: v[18751] - v[33091]]


def residuals(v, bad):
    return [fw.evalpoly(L.polys[a], v) % P for a in bad] + [f(v) % P for f in HOLD]


def solve_round(v, verbose=True):
    bad = fw.bad_checks(v)
    if not bad:
        return v, True
    r = residuals(v, bad)
    n = len(r)
    # scan controls
    cols = []
    ctrl = []
    t0 = time.time()
    for c in NONBIT:
        old = v[c]
        v[c] = old + 1
        fw.forward(v)
        r1 = residuals(v, bad)
        v[c] = old
        fw.forward(v)
        col = [(r1[i] - r[i]) % P for i in range(n)]
        if any(col):
            cols.append(col)
            ctrl.append(c)
    if verbose:
        print(f"    scan: {len(ctrl)} live controls of {len(NONBIT)} ({time.time()-t0:.0f}s)", flush=True)
    if not ctrl:
        return v, False
    J = [[cols[j][i] for j in range(len(ctrl))] for i in range(n)]
    rhs = [(-r[i]) % P for i in range(n)]
    d = gauss_solve(J, rhs, P)
    if d is None:
        if verbose:
            print("    joint system INCONSISTENT")
        return v, False
    for j, c in enumerate(ctrl):
        v[c] = (v[c] + d[j]) % P
    fw.forward(v)
    return v, True


if __name__ == '__main__':
    v, res = tri7.run(0)
    print("start:", fw.bad_checks(v), "modp:", res, flush=True)
    best = (len(L.failing_eqs(L.all_atom_values(v))), [x for x in v])
    for rnd in range(8):
        v, ok = solve_round(v)
        bad = fw.bad_checks(v)
        f = L.failing_eqs(L.all_atom_values(v))
        print(f"  round{rnd}: ok={ok} bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad[:12]}", flush=True)
        if len(f) < best[0]:
            best = (len(f), [x for x in v])
        if not bad:
            break
        if not ok:
            break
        tri7.close_all(v, set())
        bad = fw.bad_checks(v)
        f = L.failing_eqs(L.all_atom_values(v))
        print(f"    after close: bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)}", flush=True)
        if len(f) < best[0]:
            best = (len(f), [x for x in v])
        if not bad:
            break
    print(f"BEST failing={best[0]} score={L.NEQ-best[0]}")
    json.dump({str(i): best[1][i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'joint_best.json'), 'w'))
