#!/usr/bin/env python3
"""Sharp diagnostic of the tangent obstruction.

Builds the exact GF(p) Jacobian (forward AD) of ALL constraint atoms w.r.t. every
free variable in the backward cone of the violated constraints, then answers:

  Q1  which rows are inconsistent in  J d = -r ?
  Q2  rank(J_sat)  vs  rank(J_sat + one violated row)  -- is the violated residual
      in the row space of the satisfied constraints (i.e. frozen), or independent?
  Q3  restricted to ker(J_sat) (moves that disturb nothing), what is the gradient of
      each violated constraint?  Zero => locally frozen to first order.

Saves the Jacobian so later analyses need not recompute it.
"""
import sys, os, time, pickle
from collections import deque
import jengine as E, jman as J, jmodp as MP, jsolve2 as S
import jnewton2 as N

P = MP.P
HERE = os.path.dirname(os.path.abspath(__file__))


def build(b1, b2):
    cache = os.path.join(HERE, f'jjac_{b1}{b2}.pkl')
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    val, bad = S.branch(b1, b2)
    r0 = MP.residues(val)
    bad = sorted(i for i, x in r0.items() if x)
    lv, cone = N.cone_leaves(bad)
    knobs = sorted(lv)
    data = {}
    t0 = time.time()
    for jx, z in enumerate(knobs):
        for i, v in N.dcolumn(val, z).items():
            data[(i, jx)] = v
        if jx % 100 == 0:
            print(f"   col {jx}/{len(knobs)} {time.time()-t0:.0f}s", flush=True)
    obj = {'val': val, 'bad': bad, 'knobs': knobs, 'data': data, 'r0': r0}
    pickle.dump(obj, open(cache, 'wb'))
    return obj


def rowmap_of(data):
    rm = {}
    for (i, jx), v in data.items():
        rm.setdefault(i, {})[jx] = v
    return rm


def echelon(rows):
    """rows: list of (dict col->val, rhs).  Returns (piv, inconsistent_indices)."""
    piv = {}
    inc = []
    for idx, (row0, b0) in enumerate(rows):
        row = dict(row0); b = b0 % P
        while row:
            c = min(row)
            if c in piv:
                prow, pb = piv[c]
                f = row[c]
                for cc, vv in prow.items():
                    nv = (row.get(cc, 0) - f * vv) % P
                    if nv: row[cc] = nv
                    else: row.pop(cc, None)
                b = (b - f * pb) % P
            else:
                inv = pow(row[c], P - 2, P)
                piv[c] = ({cc: vv * inv % P for cc, vv in row.items()}, b * inv % P)
                row = None
                break
        if row is not None and not row and b:
            inc.append(idx)
    return piv, inc


def reduce_row(piv, row, b):
    row = dict(row); b = b % P
    while row:
        c = min(row)
        if c not in piv:
            return row, b, False          # independent
        prow, pb = piv[c]
        f = row[c]
        for cc, vv in prow.items():
            nv = (row.get(cc, 0) - f * vv) % P
            if nv: row[cc] = nv
            else: row.pop(cc, None)
        b = (b - f * pb) % P
    return row, b, True                    # dependent (row in span)


if __name__ == '__main__':
    b1, b2 = int(sys.argv[1]), int(sys.argv[2])
    obj = build(b1, b2)
    val, bad, knobs, data, r0 = obj['val'], obj['bad'], obj['knobs'], obj['data'], obj['r0']
    rm = rowmap_of(data)
    rws = sorted(rm)
    print(f"\nbranch ({b1},{b2})  violated {bad}")
    print(f"Jacobian {len(rws)} x {len(knobs)}  nnz {len(data)}")
    for i in bad:
        print(f"  a{i}: gradient support {len(rm.get(i,{}))} knobs")

    # Q1: full system consistency, naming the inconsistent rows
    allrows = [(rm[i], (-r0[i]) % P) for i in rws]
    piv, inc = echelon(allrows)
    print(f"\nQ1 full system: rank {len(piv)}, inconsistent rows {len(inc)} -> "
          f"{[rws[k] for k in inc]}")

    # Q2/Q3: echelon of the SATISFIED rows only, then test each violated row
    sat = [i for i in rws if i not in bad]
    pivS, _ = echelon([(rm[i], 0) for i in sat])
    print(f"\nQ2 rank of satisfied-only Jacobian: {len(pivS)} (of {len(knobs)} cols)"
          f"  => kernel dimension {len(knobs)-len(pivS)}")
    for i in bad:
        row, b, dep = reduce_row(pivS, rm[i], (-r0[i]) % P)
        if dep:
            print(f"  a{i}: gradient lies ENTIRELY in the span of the satisfied rows "
                  f"=> FROZEN on ker(J_sat); residual {'CANNOT' if b else 'can'} be moved"
                  f"  (reduced rhs {'nonzero' if b else 'zero'})")
        else:
            print(f"  a{i}: gradient has a component OUTSIDE the satisfied span "
                  f"=> movable at first order ({len(row)} independent cols)")

    # Q3b: joint test -- can both be fixed simultaneously?
    rows2 = [(rm[i], 0) for i in sat] + [(rm[i], (-r0[i]) % P) for i in bad]
    piv2, inc2 = echelon(rows2)
    print(f"\nQ3 satisfied(homogeneous) + violated(inhomogeneous): rank {len(piv2)}, "
          f"inconsistent {len(inc2)}")
