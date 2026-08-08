#!/usr/bin/env python3
"""ATTACK 1d (optimized): per-digit-vector FEASIBILITY of a zero-energy completion.

Faithful iff  (a zero-energy completion exists)  <=>  (dg is a true solution).
Early-exits on the first completion (feasibility only). Mersenne p=7 minimizes
the {r,r+p} representative freedom, so the search is small.  Directly tests the
degenerate x1==x2 path: colliding digit-vecs must be INFEASIBLE.
"""
import os, sys, itertools, time
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.join(HERE, '..', '..', 'squeeze')
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
sys.path.insert(0, os.path.abspath(SQ))
from ecsmall import curve, find
from ladder2 import build_win2


class Found(Exception): pass


def feasible_completion(Q, fixed):
    eqs = Q.squares
    n = Q.n
    occ = [[] for _ in range(n)]
    for e, (lin, k) in enumerate(eqs):
        for v in lin: occ[v].append(e)
    gates = [(z, i, j) for (i, j), z in Q.andcache.items()]
    gate_by_var = defaultdict(list)
    for g in gates:
        for v in g: gate_by_var[v].append(g)
    order_by_var = defaultdict(list)
    for (u, w) in Q.orders:
        order_by_var[u].append((u, w)); order_by_var[w].append((u, w))
    val = [None] * n
    for v, b in fixed.items(): val[v] = b

    def ok_local(v):
        for e in occ[v]:
            lin, k = eqs[e]
            lo = hi = k
            for vv, c in lin.items():
                if val[vv] is None: lo += min(0, c); hi += max(0, c)
                else: lo += c * val[vv]; hi += c * val[vv]
            if lo > 0 or hi < 0: return False
        for (z, i, j) in gate_by_var[v]:
            a, bb, cc = val[i], val[j], val[z]
            if a is not None and bb is not None and cc is not None and cc != a * bb: return False
            if cc == 1 and (a == 0 or bb == 0): return False
        for (u, w) in order_by_var[v]:
            if val[u] == 0 and val[w] == 1: return False
        return True

    def rec(v):
        if v == n:
            raise Found()
        if val[v] is not None:
            if ok_local(v): rec(v + 1)
            return
        for b in (0, 1):
            val[v] = b
            if ok_local(v): rec(v + 1)
        val[v] = None

    try:
        rec(0)
    except Found:
        return True
    return False


def run(p, B, m, w, mode='wallace', mult='schoolbook', leaf=3, cap_seconds=90):
    add, mul = curve(p, B)
    G, order = find(p, B)
    assert order >= (1 << m)
    M = (m + w - 1) // w
    D = 1 << w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(D)] for j in range(M)]

    def dg_point(dg):
        S = None
        for j in range(M): S = add(S, table[j][dg[j]])
        return S

    def chain_ok(dg):
        S = table[0][dg[0]]
        for j in range(1, M):
            Qp = table[j][dg[j]]
            if S is None or Qp is None or S[0] == Qp[0]: return False
            S = add(S, Qp)
        return S is not None

    kdg = next(list(dg) for dg in itertools.product(range(D), repeat=M) if chain_ok(dg))
    Tp = dg_point(kdg)
    Q, U = build_win2(p, table, Tp, w, mode=mode, mult=mult, leaf=leaf, red='naf')
    Q.finalize()
    nsol = sum(1 for dg in itertools.product(range(D), repeat=M) if dg_point(list(dg)) == Tp)
    print(f"  p={p} B={B} m={m} w={w} M={M} {mode}: vars={Q.n} sols={nsol}")
    breaks = []
    t0 = time.time()
    for dg in itertools.product(range(D), repeat=M):
        if time.time() - t0 > cap_seconds:
            print("    (time cap)"); break
        fixed = {}
        for j in range(M):
            for t in range(D): fixed[U[j][t]] = 1 if t == dg[j] else 0
        feas = feasible_completion(Q, dict(fixed))
        is_sol = dg_point(list(dg)) == Tp
        collide = not chain_ok(dg)
        st = 'SOL' if is_sol else ('degen' if collide else 'nonsol')
        flag = ''
        if feas != is_sol:
            flag = '  *** BREAK ***'; breaks.append((dg, feas, is_sol, collide))
        print(f"    dg={dg} {st:6s} feasible={feas}{flag}")
    return breaks


if __name__ == '__main__':
    sys.setrecursionlimit(1000000)
    allb = []
    for c in [(7,3,2,1,'wallace'),(7,3,2,1,'binary'),(7,1,2,1,'wallace'),(7,3,3,1,'wallace')]:
        try: allb += run(*c)
        except Exception as e:
            import traceback; print(f"  case {c} raised: {e}"); traceback.print_exc()
    print("\nTOTAL BREAKS:", len(allb), allb)
