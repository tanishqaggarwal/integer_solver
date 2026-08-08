#!/usr/bin/env python3
"""run_gadgets.py -- degree-2 SDP exactness of the encoder's arithmetic gadgets.

Builds each gadget's *actual* penalty from the read-only encoder (qubo.QB,
ladder.not_equal), brute-forces its ground energy + all ground states, then
solves the Shor SDP and reports the additive integrality gap and the rank of
the optimal moment matrix.

gap == 0  <=>  degree-2 SOS certifies the gadget (convexly exact).
"""
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import numpy as np
from qubo import QB
from ladder import not_equal
from sdp import qubo_to_C, sdp_min, rank_of, brute_min


def report(name, Q, n, extra=""):
    Emin, states = brute_min(Q, n)
    C = qubo_to_C(Q, n)
    val, X, V = sdp_min(C, restarts=32, seed=1)
    gap = Emin - val
    rk, eig = rank_of(X)
    # which variables are constant across ALL ground states (soundness gate)?
    arr = np.array(states)
    pinned = [i for i in range(n) if len(set(arr[:, i])) == 1]
    print(f"=== {name} {extra}")
    print(f"    vars={n}  ground_E={Emin}  #ground_states={len(states)}")
    print(f"    SDP_opt={val:+.6f}  additive_gap={gap:+.6f}  "
          f"{'TIGHT (deg-2 SOS exact)' if abs(gap) < 1e-4 else 'LOOSE'}")
    print(f"    rank(X*)={rk}  top_eigs={np.round(eig[:4],4)}")
    print(f"    globally-pinned vars across all ground states: "
          f"{len(pinned)} / {n}  -> {pinned if pinned else 'none'}")
    return dict(name=name, n=n, Emin=Emin, ngs=len(states), sdp=val, gap=gap,
                rank=rk, pinned=pinned)


def gadget_and():
    qb = QB()
    a = qb.new('a', 'in'); b = qb.new('b', 'in')
    qb.AND(a, b)                     # creates z, penalty ab-2az-2bz+3z
    qb.finalize()
    return qb.Q, qb.n


def gadget_fa():
    """3:2 compressor / full adder: exactly the _wallace fa penalty
       (a+b+c - s - 2d)^2  with carry-in const nc=0."""
    qb = QB()
    a = qb.new('a', 'in'); b = qb.new('b', 'in'); c = qb.new('c', 'in')
    s = qb.new('s', 'sum'); d = qb.new('d', 'carry')
    qb.add_square({a: 1, b: 1, c: 1, s: -1, d: -2}, 0)
    qb.finalize()
    return qb.Q, qb.n


def gadget_fa_nc(nc):
    """full adder with a carry-in constant nc in {1,2} (as _wallace emits when
    some inputs are the ONE literal)."""
    qb = QB()
    a = qb.new('a', 'in'); b = qb.new('b', 'in')
    s = qb.new('s', 'sum'); d = qb.new('d', 'carry')
    qb.add_square({a: 1, b: 1, s: -1, d: -2}, nc)
    qb.finalize()
    return qb.Q, qb.n


def gadget_neq(width, cst):
    qb = QB()
    word = [qb.new(f'w{i}', 'word') for i in range(width)]
    not_equal(None, qb, word, 'w', cst, f'ne{cst}')
    qb.finalize()
    return qb.Q, qb.n


if __name__ == "__main__":
    rows = []
    rows.append(report("AND  z=a*b", *gadget_and()))
    rows.append(report("FA/3:2 compressor (a+b+c-s-2d)^2", *gadget_fa()))
    rows.append(report("FA carry-in nc=1", *gadget_fa_nc(1)))
    rows.append(report("FA carry-in nc=2", *gadget_fa_nc(2)))
    for width, cst in ((3, 0), (3, 5), (4, 9)):
        rows.append(report(f"not_equal (w={width}, c={cst})", *gadget_neq(width, cst)))
    print("\n--- summary ---")
    for r in rows:
        print(f"{r['name']:38s} gap={r['gap']:+.4f} rank={r['rank']} "
              f"ngs={r['ngs']} pinned={len(r['pinned'])}")
