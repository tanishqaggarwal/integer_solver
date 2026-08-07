#!/usr/bin/env python3
"""Two exact lower bounds on failures for a support T, swept over reachable supports.

  LB1 = #{equations meeting T in exactly one atom}        (they can never cancel)
  LB2 = |R| - maxzero,  maxzero <= |T| - 1 generically, so LB2 >= |R| - |T| + 1
        (zeroing a row costs a degree of freedom; a != 0 is required)

failures >= max(LB1, |R| - |T| + 1).  Anything with that >= 7 cannot beat 39026.
Sweeps T = T0 union occ(v) over every variable v, and T0 union occ(v) union occ(w)
over the most promising pairs.
"""
import os, pickle, sys, itertools
from collections import defaultdict, Counter
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
occ = CL.occ


def bounds(T):
    T = set(T)
    if not T:
        return 0, 0, 0
    R = {i for j in T for i in atom2eq[j]}
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
    return alone, len(R), max(alone, len(R) - len(T) + 1)


if __name__ == '__main__':
    val = E.load(CL.DEL)
    DELT = CL.all_nonzero(val)
    print("deliverable:", bounds(DELT), "(alone, |R|, LB)  actual failures 7")

    bases = {
        'deliverable': DELT,
        'branch(1,1) residual': [20407, 20409, 31575],
        'branch(1,0) pins': [731, 31571],
        'on-manifold 2-atom': [8583, 30271],
        'empty': [],
    }
    allres = {}
    for lab, T0 in bases.items():
        res = []
        for v in range(E.NV):
            T = set(T0) | set(occ[v])
            if not T:
                continue
            a, r, lb = bounds(T)
            res.append((lb, a, r, len(T), v))
        res.sort()
        allres[lab] = res
        print(f"\n=== base {lab} (|T0|={len(T0)})")
        print("   best (LB, alone, |R|, |T|, var):")
        for x in res[:10]:
            print("     ", x)
        print("   #vars with LB <= 6:", sum(1 for x in res if x[0] <= 6))

    # pairs, restricted to the most promising single moves
    print("\n=== pair sweep on the deliverable base (top 60 singles) ===")
    T0 = DELT
    top = [x[4] for x in allres['deliverable'][:60]]
    best = []
    for v, w in itertools.combinations(top, 2):
        T = set(T0) | set(occ[v]) | set(occ[w])
        a, r, lb = bounds(T)
        best.append((lb, a, r, len(T), v, w))
    best.sort()
    for x in best[:12]:
        print("   ", x)

    print("\n=== pair sweep from EMPTY base (top 60 singles) ===")
    top = [x[4] for x in allres['empty'][:60]]
    best = []
    for v, w in itertools.combinations(top, 2):
        T = set(occ[v]) | set(occ[w])
        a, r, lb = bounds(T)
        best.append((lb, a, r, len(T), v, w))
    best.sort()
    for x in best[:12]:
        print("   ", x)
    pickle.dump(allres, open(os.path.join(HERE, 'jscan2.pkl'), 'wb'))
