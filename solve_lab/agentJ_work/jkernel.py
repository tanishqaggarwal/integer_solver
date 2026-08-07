#!/usr/bin/env python3
"""Beyond the tangent space.

The tangent test says: on ker(J_sat) the violated residuals have zero gradient.
But the map is polynomial (degree up to 3 in places), so a LARGE move along the
kernel can still change things at second order.  This script:
  1. computes an explicit basis of ker(J_sat) over GF(p),
  2. takes large random combinations of it, applies them EXACTLY, re-propagates,
     and reports what actually happens to the violated residues and to the rest,
  3. separately reports how much of the kernel is supported on boolean-constrained
     knobs (where a tangent test is structurally blind).
"""
import sys, os, pickle, random, time
from collections import defaultdict
import jengine as E, jman as J, jmodp as MP
import jdiag as D

P = MP.P
HERE = os.path.dirname(os.path.abspath(__file__))


def kernel_basis(rows, ncols):
    """rows: list of dict col->val (homogeneous).  Returns list of dict col->val."""
    piv = {}
    for row0 in rows:
        row = dict(row0)
        while row:
            c = min(row)
            if c in piv:
                prow = piv[c]
                f = row[c]
                for cc, vv in prow.items():
                    nv = (row.get(cc, 0) - f * vv) % P
                    if nv: row[cc] = nv
                    else: row.pop(cc, None)
            else:
                inv = pow(row[c], P - 2, P)
                piv[c] = {cc: vv * inv % P for cc, vv in row.items()}
                row = None
                break
    # back-substitute to reduced row echelon
    for c in sorted(piv):
        prow = piv[c]
        for cc in sorted([x for x in prow if x != c and x in piv]):
            f = prow.get(cc, 0)
            if not f: continue
            other = piv[cc]
            for k2, v2 in other.items():
                nv = (prow.get(k2, 0) - f * v2) % P
                if nv: prow[k2] = nv
                else: prow.pop(k2, None)
    freecols = [c for c in range(ncols) if c not in piv]
    basis = []
    for fc in freecols:
        vec = {fc: 1}
        for c, prow in piv.items():
            v = prow.get(fc, 0)
            if v:
                vec[c] = (-v) % P
        basis.append(vec)
    return basis, freecols


if __name__ == '__main__':
    b1, b2 = int(sys.argv[1]), int(sys.argv[2])
    obj = D.build(b1, b2)
    val, bad, knobs, data, r0 = obj['val'], obj['bad'], obj['knobs'], obj['data'], obj['r0']
    rm = D.rowmap_of(data)
    rws = sorted(rm)
    sat = [i for i in rws if i not in bad]
    print(f"branch ({b1},{b2}) violated {bad}; rows {len(rws)} knobs {len(knobs)}")

    basis, freecols = kernel_basis([rm[i] for i in sat], len(knobs))
    print(f"ker(J_sat) dimension {len(basis)}")

    # which knobs carry a boolean constraint?
    boolknob = set()
    for jx, z in enumerate(knobs):
        for i in MP.CONS:
            if E.varsof[i] == {z} and any(len(k) == 2 for k in E.polys[i]):
                boolknob.add(jx); break
    print(f"knobs carrying a boolean constraint: {len(boolknob)} of {len(knobs)}")

    onlybool = 0
    for vec in basis:
        sup = set(vec)
        if sup <= boolknob:
            onlybool += 1
    print(f"kernel basis vectors supported ONLY on boolean knobs: {onlybool}/{len(basis)}")

    # large nonlinear probes along the kernel
    random.seed(99)
    print("\nlarge random kernel moves (exact re-propagation):")
    for trial in range(12):
        v2 = list(val)
        for vec in basis:
            c = random.randrange(P)
            if c == 0: continue
            for jx, coef in vec.items():
                v2[knobs[jx]] = (v2[knobs[jx]] + c * coef) % P
        MP.fwd_modp(v2)
        r = MP.residues(v2)
        nb = sorted(i for i, x in r.items() if x)
        same = all(r[i] == r0[i] for i in bad)
        print(f"  trial {trial}: violated {len(nb)}  residues-of-violated-unchanged={same}"
              f"  newly_broken={len([i for i in nb if i not in bad])}")

    # single-basis-vector moves, large coefficient
    print("\nsingle kernel direction, large coefficient:")
    hits = 0
    for k, vec in enumerate(basis):
        v2 = list(val)
        c = random.randrange(P)
        for jx, coef in vec.items():
            v2[knobs[jx]] = (v2[knobs[jx]] + c * coef) % P
        MP.fwd_modp(v2)
        r = MP.residues(v2)
        nb = [i for i in MP.CONS if r[i]]
        moved = [i for i in bad if r[i] != r0[i]]
        if moved:
            hits += 1
            print(f"  basis {k}: MOVES {moved}; total violated {len(nb)}")
    print(f"kernel directions that move a violated residue at all: {hits}/{len(basis)}")
