#!/usr/bin/env python3
"""Can free-knob atoms cancel a reachable residual?

Base residual T0 (values fixed, constrained) + u free-knob atoms (values ANY integer,
zero side effects).  Rows R = eqs(T0 union K).  Over the free-knob coordinates the
lattice is all of Z^u, so integrality is free: any rational kernel vector scales to an
integer one.  Hence

    min failures = |R| - max #rows simultaneously zeroable

and with the T0 coordinates pinned we can zero a set S of rows iff the linear system
  sum_{alpha in K} c_{i,alpha} y_alpha = -sum_{j in T0} c_{ij} a_j   (i in S)
has a rational solution, i.e. iff rank(C_S,K) == rank([C_S,K | -C_S,T0 a]).

This script measures the overlap structure that decides it.
"""
import os, pickle, sys, itertools
from fractions import Fraction
from collections import defaultdict, Counter
import jengine as E
import jcluster as CL
import jexpand as X

atom2eq, eqrow, occ = X.atom2eq, X.eqrow, X.occ
KNOBS = set(X.knobs) if hasattr(X, 'knobs') else None
if KNOBS is None:
    KNOBS = {i for i in range(len(E.polys)) if X.private_var(i)}


def rank_aug(rows_S, cols, rhs):
    """rational rank of C and of [C|rhs]."""
    A = [[Fraction(eqrow[i].get(c, 0)) for c in cols] + [Fraction(rhs[k])]
         for k, i in enumerate(rows_S)]
    m, n = len(A), len(cols)
    r = 0
    for c in range(n):
        p = None
        for i in range(r, m):
            if A[i][c]:
                p = i; break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [A[i][k] - f * A[r][k] for k in range(n + 1)]
        r += 1
    # rank of augmented
    ra = r
    for i in range(r, m):
        if A[i][n]:
            ra += 1
            break
    return r, ra


def analyse(T0, a0, label):
    R0 = sorted({i for j in T0 for i in atom2eq[j]})
    R0set = set(R0)
    print(f"\n=== {label}: |T0|={len(T0)} |R0|={len(R0)}")
    # free knobs ranked by how many NEW rows they bring
    cands = []
    for al in KNOBS:
        new = [i for i in atom2eq[al] if i not in R0set]
        inside = len(atom2eq[al]) - len(new)
        if inside:
            cands.append((len(new), inside, al))
    cands.sort()
    print(f"  free-knob atoms sharing >=1 equation with R0: {len(cands)}")
    for c in cands[:15]:
        print(f"     new_rows={c[0]:3d} shared={c[1]:3d} atom a{c[2]}")
    if not cands:
        print("  NONE -> free knobs cannot touch this residual at all")
        return
    # greedy: take knobs that add zero new rows
    zero_new = [c[2] for c in cands if c[0] == 0]
    print(f"  knobs adding ZERO new rows: {len(zero_new)} -> {zero_new[:20]}")
    if zero_new:
        rhs = [-sum(eqrow[i].get(j, 0) * a0[k] for k, j in enumerate(T0)) for i in R0]
        best = 0
        for u in range(1, min(len(zero_new), 12) + 1):
            K = zero_new[:u]
            # max rows zeroable = largest S with rank consistency; test all rows at once
            r, ra = rank_aug(R0, K, rhs)
            print(f"    u={u}: rank {r}, augmented {ra}  -> all {len(R0)} rows "
                  f"{'CONSISTENT' if r == ra else 'inconsistent'}")


if __name__ == '__main__':
    val = E.load(CL.DEL)
    T0 = CL.all_nonzero(val)
    a0 = [CL.atomval(i, val) for i in T0]
    analyse(T0, a0, "deliverable support")

    for T0x, lab in [([8583, 30271], 'on-manifold 2-atom'),
                     ([731, 31571], 'branch (1,0) pins'),
                     ([20407, 20409, 31575], 'branch (1,1) residual')]:
        analyse(T0x, [1] * len(T0x), lab)
