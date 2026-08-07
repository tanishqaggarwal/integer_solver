#!/usr/bin/env python3
"""Atom-level minimum-distance analysis.

SCOPE FIRST.  score = 39033 - ||M a||_0 with a = atom(x).  One might hope to bound the
score by relaxing to all nonzero a in Z^n.  That relaxation is VACUOUS, unconditionally:

  * if rank_Q(M) < n, some nonzero a has M a = 0            -> relaxed min weight 0
  * if rank_Q(M) = n, then adj(M) e_1 is a nonzero integer vector and
    M (adj(M) e_1) = det(M) e_1                             -> relaxed min weight 1

so min over nonzero a in Z^n of ||M a||_0 <= 1 either way, and no rank computation can
change that.  Any bound above 1 MUST constrain the achievable set.  The right currency
is therefore the support, which is what the un-cancellable-row bound uses:

  ||M a||_0 >= alone(T) := #{equations meeting T = supp(a) in exactly one atom}   [exact]

and, when every (|T|-1)-subset of R(T) has full rank, also >= |R(T)| - |T| + 1
(heuristic, flagged as such -- it is the generic-position count, not a theorem).

This script searches for supports that could host a weight-<7 codeword, by beam search
on |R(T)| - |T| (adding an atom costs +1 to |T| and +|eqs(j)\\R| to |R|, so an atom
bringing 0 new rows strictly improves the block).
"""
import os, pickle, sys, random, heapq
from collections import defaultdict
import jengine as E
import jcluster as CL

HERE = os.path.dirname(os.path.abspath(__file__))
atom2eq = defaultdict(list)
eqrow = []
for e in CL.eqs:
    row = {}
    for c, j in e['terms']:
        row[j] = row.get(j, 0) + c
    row = {j: c for j, c in row.items() if c}
    eqrow.append(row)
    for j in row:
        atom2eq[j].append(e['i'])
A2E = {j: frozenset(v) for j, v in atom2eq.items()}
NA = len(eqrow)


def stats(T):
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
    return alone, len(R), len(T)


def candidates(T, R):
    c = set()
    for i in R:
        c |= set(eqrow[i])
    return c - T


def beam(seed, width=40, steps=40):
    """beam search minimising |R| - |T|; returns best (|R|-|T|, alone, |T|, T)."""
    T0 = frozenset(seed)
    R0 = set()
    for j in T0:
        R0 |= A2E[j]
    best = (len(R0) - len(T0), None, T0)
    beam_ = [(len(R0) - len(T0), T0, frozenset(R0))]
    seen = {T0}
    for _ in range(steps):
        nxt = []
        for score, T, R in beam_:
            for j in candidates(set(T), R):
                nR = R | A2E[j]
                nT = T | {j}
                if nT in seen:
                    continue
                seen.add(nT)
                nxt.append((len(nR) - len(nT), nT, frozenset(nR)))
        if not nxt:
            break
        nxt.sort(key=lambda t: (t[0], len(t[1])))
        beam_ = nxt[:width]
        if beam_[0][0] < best[0]:
            best = (beam_[0][0], None, beam_[0][1])
    a, r, t = stats(best[2])
    return best[0], a, r, t, best[2]


if __name__ == '__main__':
    val = E.load(CL.DEL)
    DELT = CL.all_nonzero(val)
    a, r, t = stats(DELT)
    print(f"deliverable cluster: |T|={t} |R|={r} alone={a} |R|-|T|={r-t} "
          f"(actual weight 7)")

    print("\nbeam search from the deliverable cluster:")
    d, a, r, t, T = beam(DELT, width=60, steps=40)
    print(f"  best |R|-|T| = {d}  (|T|={t} |R|={r} alone={a})")
    print(f"  support size {len(T)}")

    print("\nbeam search from each atom of the cluster:")
    for j in DELT:
        d, a, r, t, T = beam([j], width=30, steps=30)
        print(f"  seed a{j}: |R|-|T|={d} |T|={t} |R|={r} alone={a}")

    random.seed(5)
    print("\nbeam search from 60 random atom seeds:")
    best = None
    for j in random.sample(range(NA), 60):
        d, a, r, t, T = beam([j], width=20, steps=25)
        if best is None or (d, a) < (best[0], best[1]):
            best = (d, a, r, t, j)
    print("  best over random seeds:", best)
