#!/usr/bin/env python3
"""Can a cluster of confined pins beat weight 7?

Break a set T of pins.  Their pin atoms are the support; their variables are knobs.
Three degrees of freedom are consumed zeroing the residual constraints, so

    weight >= max( alone(T),  |R(T)| - (|T| - 3) )

Both bounds are computed exactly (alone(T) is the un-cancellable-row bound).
Beam-search the 83 confined pins for the cluster minimising that.
"""
import os, sys, pickle, itertools
from collections import defaultdict
import jengine as E
import jdist as DI

HERE = os.path.dirname(os.path.abspath(__file__))
A2E, eqrow = DI.A2E, DI.eqrow

cand = pickle.load(open(os.path.join(HERE, 'jpins2.pkl'), 'rb'))
cand.sort()
PINS = [(c, v, p, m) for (c, v, p, m) in cand]
print(f"confined pins: {len(PINS)}")
pinatoms = [p for (c, v, p, m) in PINS]


def bounds(Tatoms):
    T = set(Tatoms)
    R = set()
    for j in T:
        R |= A2E[j]
    alone = 0
    for i in R:
        n = 0
        for j in eqrow[i]:
            if j in T:
                n += 1
                if n > 1:
                    break
        if n == 1:
            alone += 1
    dof = max(0, len(T) - 3)
    return alone, len(R), max(alone, len(R) - dof)


if __name__ == '__main__':
    print("\nsingle pins (alone, |R|, bound):")
    for c, v, p, m in PINS[:6]:
        print(f"   a{p} cost {c}: {bounds([p])}")

    # beam search over subsets of the confined pin atoms
    best = None
    beam = [frozenset([p]) for p in pinatoms]
    seen = set(beam)
    for step in range(12):
        scored = []
        for Tt in beam:
            a, r, b = bounds(Tt)
            scored.append((b, r - max(0, len(Tt) - 3), len(Tt), Tt))
        scored.sort(key=lambda t: (t[0], t[1]))
        if best is None or scored[0][0] < best[0]:
            best = (scored[0][0], scored[0][3])
        keep = [t[3] for t in scored[:120]]
        nxt = []
        for Tt in keep:
            for p in pinatoms:
                if p in Tt:
                    continue
                nT = Tt | {p}
                if nT in seen:
                    continue
                seen.add(nT); nxt.append(nT)
        if not nxt:
            break
        beam = nxt
        a, r, b = bounds(scored[0][3])
        print(f"  step {step}: best bound {scored[0][0]} (|T|={len(scored[0][3])} "
              f"|R|={r} alone={a})")
    a, r, b = bounds(best[1])
    print(f"\nBEST over confined-pin clusters: bound {best[0]}  "
          f"|T|={len(best[1])} |R|={r} alone={a}")
    print("  support atoms:", sorted(best[1]))
    print(f"  => best achievable score by this route <= {39033 - best[0]}")
    print(f"  (deliverable = 39026, i.e. weight 7)")

    # explicit families that look tight
    fams = [[4795, 4796, 4797, 4798, 4799], [12151, 12152, 12153, 12154],
            [29698, 29699, 29700, 29701, 29702], [35445, 35446, 35447, 35448, 35449],
            [638, 639, 640, 641], [39032]]
    print("\nnamed families:")
    for f in fams:
        f2 = [x for x in f if x in set(pinatoms)]
        if f2:
            print(f"   {f2}: {bounds(f2)}")
    print("\nfamily + a39032:")
    for f in fams[:-1]:
        f2 = [x for x in f if x in set(pinatoms)] + [39032]
        print(f"   {f2}: {bounds(f2)}")
