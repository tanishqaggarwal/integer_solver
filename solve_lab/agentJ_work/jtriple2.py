#!/usr/bin/env python3
"""Corrected accounting for the confined-pin route.

Let S be the subset of the 3 residual constraints LEFT violated (so they join the
support), and let K be the set of pins broken (|K| = k).  Zeroing the other 3 - |S|
residual constraints consumes 3 - |S| degrees of freedom, which requires
k >= 3 - |S|.  Support T = K u S, and

    weight >= max( alone(T),  |R(T)| - (k - (3 - |S|)) ).

The previous script credited |T| < 3 with dof it does not have; this one does not.
"""
import os, sys, pickle, itertools
from collections import defaultdict
import jengine as E
import jdist as DI

HERE = os.path.dirname(os.path.abspath(__file__))
A2E, eqrow = DI.A2E, DI.eqrow
RES = [20407, 20409, 31575]

cand = pickle.load(open(os.path.join(HERE, 'jpins2.pkl'), 'rb'))
cand.sort()
pinatoms = [p for (c, v, p, m) in cand]
# which residual constraints each pin can move (first-order census, exact-verified
# for the 40 cheapest)
moves = {p: set(m) for (c, v, p, m) in cand}
print(f"confined pins: {len(pinatoms)}")


def alone_R(T):
    T = set(T)
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
    return alone, len(R)


def bound(K, S):
    need = 3 - len(S)
    k = len(K)
    if k < need:
        return None
    # the pins must jointly be able to move the constraints being zeroed
    cover = set()
    for p in K:
        cover |= moves[p]
    if not set(RES) - set(S) <= cover:
        return None
    T = set(K) | set(S)
    alone, R = alone_R(T)
    dof = k - need
    return max(alone, R - dof), R, alone, len(T)


if __name__ == '__main__':
    best = None
    # exhaustive over S, and beam over pin sets for each S
    for r in range(4):
        for S in itertools.combinations(RES, r):
            need = 3 - len(S)
            beam = [frozenset([p]) for p in pinatoms]
            seen = set(beam)
            for step in range(10):
                scored = []
                for K in beam:
                    b = bound(K, S)
                    if b:
                        scored.append((b[0], b[1], K))
                if scored:
                    scored.sort()
                    if best is None or scored[0][0] < best[0]:
                        best = (scored[0][0], set(scored[0][2]), set(S))
                        print(f"  new best bound {best[0]}  S={S} |K|={len(best[1])}")
                keep = [K for (_, _, K) in sorted(
                    [(len(set().union(*[A2E[p] for p in K])) - len(K), 0, K)
                     for K in beam])[:80]]
                nxt = []
                for K in keep:
                    for p in pinatoms:
                        if p in K:
                            continue
                        nK = K | {p}
                        if nK in seen:
                            continue
                        seen.add(nK); nxt.append(nK)
                if not nxt:
                    break
                beam = nxt
    print(f"\nBEST over the confined-pin route: weight >= {best[0]}")
    print(f"   leave violated S = {sorted(best[2])}, break pins {sorted(best[1])}")
    print(f"   => score <= {39033 - best[0]}   (deliverable 39026)")
