#!/usr/bin/env python3
"""Enlarge the support with FREE-KNOB atoms instead of enlarging the lattice.

An atom is a free knob if it contains a variable that occurs in NO other atom and
appears linearly with coefficient +-1: setting that variable moves this atom's value
to any integer whatsoever and changes nothing else.  Such an atom contributes a fully
unconstrained coordinate (all of Z), not a coset of a q-multiple lattice.

Adding a free knob alpha to the support T costs |eqs(alpha) \\ R| new rows and buys one
unconstrained unknown.  If we can add knobs that bring few or no new rows, the number
of cancellable rows rises and the score with it:
      score = 39033 - ||M a||_0 ,  and with u unconstrained coordinates we can zero
      any set of rows whose restricted matrix has rank <= u-1 (scale to clear
      denominators -- integrality is free once a coordinate ranges over all of Z).
"""
import os, pickle, sys, itertools
from collections import defaultdict, Counter
import jengine as E
import jcluster as CL

HERE = os.path.dirname(os.path.abspath(__file__))
polys = E.polys
NA = len(polys)

# equations per atom, and rows
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


def private_var(i):
    """a variable occurring only in atom i, linearly, coefficient +-1 -> free knob."""
    p = polys[i]
    higher = set()
    for k in p:
        if len(k) >= 2:
            higher.update(k)
    for k, c in p.items():
        if len(k) == 1 and abs(c) == 1 and k[0] not in higher and len(occ[k[0]]) == 1:
            return k[0], c
    return None


if __name__ == '__main__':
    knobs = {}
    for i in range(NA):
        pv = private_var(i)
        if pv:
            knobs[i] = pv
    print(f"atoms that are FREE KNOBS (private occurrence-1 linear variable): {len(knobs)}")
    print("equations-per-knob-atom histogram:",
          Counter(len(atom2eq[i]) for i in knobs).most_common(8))

    val = E.load(CL.DEL)
    T = CL.all_nonzero(val)
    R = sorted({i for j in T for i in atom2eq[j]})
    Rset = set(R)
    print(f"\ndeliverable support {T}, rows R = {R}")

    # every atom appearing in those 12 equations
    AR = sorted({j for i in R for j in eqrow[i]})
    print(f"atoms appearing in R: {len(AR)}")
    cand = []
    for j in AR:
        if j in T:
            continue
        outside = [i for i in atom2eq[j] if i not in Rset]
        cand.append((len(outside), j, j in knobs, len(atom2eq[j])))
    cand.sort()
    print("\ncandidate atoms to add (new_rows, atom, is_free_knob, #eqs):")
    for c in cand[:40]:
        print("   ", c)
    nk = [c for c in cand if c[2]]
    print(f"\nfree-knob candidates: {len(nk)}; cheapest new-row counts: {[c[0] for c in nk[:20]]}")

    # global search: any free-knob atom whose equations are few
    best = sorted(knobs, key=lambda i: len(atom2eq[i]))[:20]
    print("\nglobally sparsest free-knob atoms (atom, #eqs):",
          [(i, len(atom2eq[i])) for i in best])
