#!/usr/bin/env python3
"""demo_win2.py -- the squeezed ladder, checked by the lab's own exhaustive test.

Identical protocol to ../demo_win.py: build the real Hamiltonian for a scaled
instance, fill every ancilla by replaying the construction, and enumerate EVERY
candidate scalar.  The encoding is faithful iff

    { k : E(k) = 0 }  is exactly  { k : k*G = T }  on the small curve.

The curves are chosen with pseudo-Mersenne p so that the NAF reduction is being
exercised, not bypassed.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ecsmall import curve, find                                  # noqa: E402
from ladder2 import build_win2                                   # noqa: E402


def make(p, B, m, w, k):
    add, mul = curve(p, B)
    G, order = find(p, B)
    assert order >= (1 << m), f"curve too small: order {order}"
    M = (m + w - 1) // w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(1 << w)]
             for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))
    Tp = add(mul(k % order, G), mul(off % order, G))
    return G, order, Tp, table, M, add, mul


def chain_ok(table, M, dg, add):
    S = table[0][dg[0]]
    for j in range(1, M):
        Q = table[j][dg[j]]
        if S is None or Q is None or S[0] == Q[0]:
            return False
        S = add(S, Q)
    return S is not None


def run(p, B, m, w, k, mode='wallace', mult='karatsuba', leaf=4, red='naf',
        naf_consts=True, exhaustive=True):
    dgs = lambda kk: [(kk >> (w * j)) % (1 << w)                 # noqa: E731
                      for j in range((m + w - 1) // w)]
    G, order, Tp, table, M, add, mul = make(p, B, m, w, k)
    if not chain_ok(table, M, dgs(k), add):
        k = next(kk for kk in range(1 << m) if chain_ok(table, M, dgs(kk), add))
        G, order, Tp, table, M, add, mul = make(p, B, m, w, k)
    t0 = time.time()
    Q, U = build_win2(p, table, Tp, w, mode=mode, mult=mult, leaf=leaf, red=red,
                      naf_consts=naf_consts)
    st = Q.stats()
    tag = (f"p={p} m={m} w={w} windows={M} {mode}/{mult}(leaf={leaf})/{red}"
           f"{'/naf-consts' if naf_consts else ''}")
    print(f"  {tag}\n    {st['vars']} vars, {st['couplers']} couplers, "
          f"K={st['max_clique']}, |J|=2^{st['dynamic_range_bits']}, "
          f"{time.time() - t0:.1f}s build")

    def energy_of(kk):
        dg = dgs(kk)
        wv0 = {f"_u{j}": dg[j] for j in range(M)}
        inp = {}
        for j in range(M):
            for t in range(1 << w):
                inp[U[j][t]] = 1 if t == dg[j] else 0
        try:
            x, _ = Q.witness(inp, wv0)
        except (AssertionError, ValueError, ZeroDivisionError):
            return None
        return Q.energy(x)

    if energy_of(k) != 0:
        for kk in range(1 << m):
            if mul(kk, G) is not None and energy_of(kk) == 0:
                k = kk
                break
    assert energy_of(k) == 0, "true k is not a zero-energy state"
    if exhaustive:
        zeros = [kk for kk in range(1 << m) if energy_of(kk) == 0]
        degen = sum(1 for kk in range(1 << m) if energy_of(kk) is None)
        T = mul(k, G)
        sol = [kk for kk in range(1 << m) if mul(kk, G) == T]
        ok = set(zeros) <= set(sol) and k in zeros
        print(f"    zero-energy k: {zeros}  true solutions: {sol}  "
              f"degenerate-chain k: {degen}/{1 << m}  -- "
              f"{'FAITHFUL' if ok else '*** NOT FAITHFUL ***'}")
        assert ok
    return st


if __name__ == '__main__':
    print("=" * 78)
    print("EXHAUSTIVE FAITHFULNESS OF THE SQUEEZED LADDER  (demo_win.py protocol)")
    print("=" * 78)
    # (p, B, m) with p pseudo-Mersenne so the NAF reduction is exercised
    CASES = [(127, 3, 4), (251, 1, 4), (1021, 3, 5)]
    for p, B, m in CASES:
        for w in (1, 2, 3):
            for mode in ('wallace', 'binary'):
                for mult in ('schoolbook', 'karatsuba'):
                    for nc in (True, False):
                        run(p, B, m, w, 9, mode=mode, mult=mult, leaf=4,
                            naf_consts=nc)
    print()
    print("all variants faithful")
