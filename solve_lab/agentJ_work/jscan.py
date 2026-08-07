#!/usr/bin/env python3
"""Exhaustive scan over reachable supports.

From a base state whose residual support is T0, perturbing a single variable v makes
exactly the atoms containing v nonzero, so T = T0 union occ(v) is reachable for every
variable v.  For every such T compute the exact support-only lower bound

    failures >= #{equations meeting T in exactly one atom}

which is a NECESSARY condition: any T whose bound exceeds 6 can never beat 39026.
Pure set arithmetic, so all 38748 variables can be swept.
"""
import os, pickle, sys
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


def bound(T):
    T = set(T)
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
    return alone, len(R)


if __name__ == '__main__':
    val = E.load(CL.DEL)
    DELT = CL.all_nonzero(val)
    bases = {
        'empty (full solution)': [],
        'deliverable': DELT,
        'on-manifold 2-atom': [8583, 30271],
        'branch(1,0) pins': [731, 31571],
        'branch(1,1) residual': [20407, 20409, 31575],
    }
    for lab, T0 in bases.items():
        b0, r0 = bound(T0) if T0 else (0, 0)
        print(f"\n=== base {lab}: T0={T0} bound={b0} |R|={r0}")
        res = []
        for v in range(E.NV):
            T = set(T0) | set(occ[v])
            if not T:
                continue
            a, r = bound(T)
            res.append((a, r, v, len(occ[v])))
        res.sort()
        print("  best single-variable perturbations (bound, |R|, var, #atoms):")
        for x in res[:12]:
            print("     ", x)
        print("  #variables with bound <= 6:", sum(1 for x in res if x[0] <= 6))
        print("  bound histogram:", Counter(x[0] for x in res).most_common(8))
        pickle.dump(res, open(os.path.join(HERE, f'jscan_{abs(hash(lab))%10000}.pkl'), 'wb'))
