#!/usr/bin/env python3
"""Decisive structural diagnostic for the 4 obstruction checks.

Questions:
 1) At all-0, which atoms violate (Z)?  (expect 1817,30378,40782,44271)
 2) Structure of the 4 checks: vars, degrees.
 3) Values of x_9770,x_3183,x_18274,x_17728 at all-0.
 4) For each control bit, single-flip: which of the 4 checks' RESIDUALS move,
    and do x_9770/x_3183 (22-side) vs x_18274/x_17728 (233-side) respond?
 5) Are the 4 check residuals LINEAR in the control bits? (test random pairs
    and a random large subset vs the single-flip linear prediction)
"""
import json, time, sys
from confluent_eval5 import build5, make_forward
from propagate import NVARS

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
CHECKS = [1817, 30378, 40782, 44271]
WATCH = [9770, 3183, 18274, 17728]

def main():
    t0 = time.time()
    A, kind, info, seq, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in set(BITS22)]
    print(f"built. control={len(control)} bits22={len(BITS22)} bits233={len(bits233)} ({time.time()-t0:.0f}s)", flush=True)

    def resid(val, a):
        s = 0
        for m, c in A[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        return s

    base = solve(list(bestval), [])
    viol = [a for a in range(len(A)) if resid(base, a) != 0]
    print(f"forward_Z([]) violated atoms ({len(viol)}): {sorted(viol)}", flush=True)
    for a in CHECKS:
        vs = sorted(set().union(*[set(m) for m in A[a]])) if A[a] else []
        degs = sorted(set(len(m) for m in A[a]))
        print(f"  atom {a}: nvars={len(vs)} degs={degs} nterms={len(A[a])} resid0={resid(base,a)}", flush=True)
    print("watch values at all-0:", flush=True)
    for w in WATCH:
        print(f"  x_{w} = {base[w]}", flush=True)
    R0 = {a: resid(base, a) for a in CHECKS}

    # single-flip responses
    print("\nsingle-flip scan over control bits...", flush=True)
    d = {a: {} for a in CHECKS}          # atom -> {bit: delta_resid}
    wmove = {w: [] for w in WATCH}        # watch var -> bits that move it
    for i, b in enumerate(control):
        val = solve(list(bestval), [b])
        for a in CHECKS:
            dr = resid(val, a) - R0[a]
            if dr != 0: d[a][b] = dr
        for w in WATCH:
            if val[w] != base[w]: wmove[w].append(b)
    for a in CHECKS:
        movers = list(d[a].keys())
        n22 = sum(1 for b in movers if b in set(BITS22))
        n233 = len(movers) - n22
        print(f"  check {a}: moved by {len(movers)} bits ({n22} of the 22, {n233} of the 233)", flush=True)
    for w in WATCH:
        n22 = sum(1 for b in wmove[w] if b in set(BITS22))
        print(f"  x_{w} moved by {len(wmove[w])} bits ({n22} of 22, {len(wmove[w])-n22} of 233)", flush=True)

    # linearity test of check residuals over the control bits
    print("\nlinearity test of check residuals...", flush=True)
    import itertools
    # a) pairs among movers of atom 1817
    tests = []
    m1817 = list(d[1817].keys())
    if len(m1817) >= 4:
        tests += [list(p) for p in itertools.combinations(m1817[:6], 2)][:5]
    # b) a big subset: every other of atom-1817 movers
    if len(m1817) >= 6:
        tests.append(m1817[::2])
    # c) mixed 22+233
    mm = [b for b in m1817 if b in set(BITS22)][:3] + [b for b in m1817 if b not in set(BITS22)][:3]
    if len(mm) >= 2: tests.append(mm)
    for S in tests:
        val = solve(list(bestval), S)
        for a in CHECKS:
            actual = resid(val, a) - R0[a]
            pred = sum(d[a].get(b, 0) for b in S)
            tag = "LIN" if actual == pred else "NONLIN"
            if a == 1817 or tag == "NONLIN":
                print(f"  S(|{len(S)}|) atom {a}: {tag}  actual-pred={actual-pred}", flush=True)
    print(f"\ndone ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
