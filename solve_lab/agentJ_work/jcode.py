#!/usr/bin/env python3
"""Off-manifold coding analysis.

Score = #{i : S_i = 0}, S_i = sum_j c_ij a_j over the atom-value vector a.
So failures = ||M a||_0 and we want the minimum-weight syndrome over REALIZABLE a.

Hard fact used throughout: if atom j is in the support of a then a_j != 0, so any
equation whose only support-atom is j necessarily fails.  That gives an exact,
support-only lower bound  failures >= #{i : |supp(a) ∩ atoms(i)| == 1}  which needs
no knowledge of the values at all.
"""
import os, pickle, sys, json
from collections import defaultdict, Counter
import jengine as E

HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
eqs = M['eqs']
NA = len(M['atoms'])

# atom -> equations, equation -> {atom: coef}
atom2eq = defaultdict(list)
eqrow = []
for e in eqs:
    row = {}
    for c, j in e['terms']:
        row[j] = row.get(j, 0) + c
    row = {j: c for j, c in row.items() if c}
    eqrow.append(row)
    for j in row:
        atom2eq[j].append(e['i'])


def analyse(T, label=''):
    T = set(T)
    R = sorted({i for j in T for i in atom2eq[j]})
    only1 = []
    multi = []
    for i in R:
        s = [j for j in eqrow[i] if j in T]
        if len(s) == 1:
            only1.append(i)
        else:
            multi.append(i)
    print(f"\n=== {label}  |T|={len(T)}")
    print(f"  equations touched |R| = {len(R)}")
    print(f"  rows with exactly ONE support atom (ALWAYS fail) = {len(only1)}")
    print(f"  rows with >=2 support atoms (cancellable)        = {len(multi)}")
    print(f"  => exact support-only lower bound on failures    = {len(only1)}")
    return R, only1, multi


def per_atom(T):
    T = set(T)
    print("  per-atom: (atom, #eqs, #eqs where it is alone in T)")
    for j in sorted(T):
        alone = sum(1 for i in atom2eq[j] if len([k for k in eqrow[i] if k in T]) == 1)
        print(f"    a{j}: {len(atom2eq[j])} eqs, alone in {alone}")


if __name__ == '__main__':
    # support of the deliverable (measured, not assumed)
    P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
    polys = P['polys']
    val = E.load(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
    av = []
    for p in polys:
        s = 0
        for k, c in p.items():
            t = c
            for j in k:
                t *= val[j]
            s += t
        av.append(s)
    DEL = [i for i, x in enumerate(av) if x]
    print("deliverable support:", DEL)
    analyse(DEL, 'DELIVERABLE (39026)')
    per_atom(DEL)

    for T, lab in [([8583, 30271], 'on-manifold CD-fixed (2 atoms)'),
                   ([731, 31571], 'branch (1,0) output pins'),
                   ([20407, 20409, 31575], 'branch (1,1) nonlinear residual'),
                   ([731, 3895, 31571], 'lifted J_b10_38998')]:
        analyse(T, lab)
        per_atom(T)

    # global: distribution of #equations per atom
    print("\nequations-per-atom histogram:",
          Counter(len(v) for v in atom2eq.values()).most_common(8))
    print("atoms in <=9 equations:", sum(1 for v in atom2eq.values() if len(v) <= 9))
