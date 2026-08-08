#!/usr/bin/env python3
"""ATTACK 1d: the DECISIVE ladder faithfulness test.

demo_win2 only replays ONE canonical witness per digit-vector, so it cannot
prove a non-solution digit-vector has NO zero-energy completion.  Here we clamp
the one-hot digit inputs to each digit-vector in turn and enumerate ALL E=0
completions of the remaining (arithmetic) variables with the validated
zero_states search.  Faithful iff:
   #completions(dg) > 0   <=>   dg is a true discrete-log solution.
A non-solution digit-vector with >=1 completion is a SPURIOUS GROUND STATE.
This also directly exercises the degenerate x1==x2 (d==0) path: any wrong dg
whose chain collides must show 0 completions (the d!=0 gadget must block it).
"""
import os, sys, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.join(HERE, '..', '..', 'squeeze')
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
sys.path.insert(0, os.path.abspath(SQ))
from ecsmall import curve, find
from ladder2 import build_win2
import verify


def zero_states_clamped(Q, fixed):
    """zero_states with some variables pre-fixed to given values (dict var->bit)."""
    eqs = Q.squares
    gates = [(z, i, j) for (i, j), z in Q.andcache.items()]
    orders = Q.orders
    n = Q.n
    occ = [[] for _ in range(n)]
    for e, (lin, k) in enumerate(eqs):
        for v in lin:
            occ[v].append(e)
    out = []
    val = [None] * n
    for v, b in fixed.items():
        val[v] = b

    def feasible(elist):
        for e in elist:
            lin, k = eqs[e]
            lo = hi = k
            for v, c in lin.items():
                if val[v] is None:
                    lo += min(0, c); hi += max(0, c)
                else:
                    lo += c * val[v]; hi += c * val[v]
            if lo > 0 or hi < 0:
                return False
        return True

    # check gates/orders consistency incrementally
    def rec(v):
        if v == n:
            out.append(list(val)); return
        if val[v] is not None:  # pre-fixed
            if feasible(occ[v]) and gates_ok(v) and orders_ok(v):
                rec(v + 1)
            return
        for b in (0, 1):
            val[v] = b
            if feasible(occ[v]) and gates_ok(v) and orders_ok(v):
                rec(v + 1)
        val[v] = None

    def gates_ok(v):
        for (z, i, j) in gates:
            a, bb, cc = val[i], val[j], val[z]
            if a is not None and bb is not None and cc is not None and cc != a * bb:
                return False
            if cc == 1 and ((a == 0) or (bb == 0)):
                return False
        return True

    def orders_ok(v):
        for (u, w) in orders:
            if val[u] == 0 and val[w] == 1:
                return False
        return True

    rec(0)
    return out


def run(p, B, m, w, mode='wallace', mult='schoolbook', leaf=3, cap_seconds=60):
    add, mul = curve(p, B)
    G, order = find(p, B)
    assert order >= (1 << m)
    M = (m + w - 1) // w
    D = 1 << w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(D)] for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))

    def dg_point(dg):
        S = None
        for j in range(M):
            S = add(S, table[j][dg[j]])
        return S

    # choose target from a non-degenerate solution scalar
    def chain_ok(dg):
        S = table[0][dg[0]]
        for j in range(1, M):
            Q = table[j][dg[j]]
            if S is None or Q is None or S[0] == Q[0]:
                return False
            S = add(S, Q)
        return S is not None
    kdg = next((list(dg) for dg in itertools.product(range(D), repeat=M) if chain_ok(dg)), None)
    assert kdg is not None
    Tp = dg_point(kdg)

    Q, U = build_win2(p, table, Tp, w, mode=mode, mult=mult, leaf=leaf, red='naf')
    Q.finalize()
    print(f"  p={p} B={B} m={m} w={w} M={M} D={D} {mode}/{mult}: vars={Q.n}, |target sols|="
          f"{sum(1 for dg in itertools.product(range(D),repeat=M) if dg_point(list(dg))==Tp)}")

    breaks = []
    t0 = time.time()
    for dg in itertools.product(range(D), repeat=M):
        if time.time() - t0 > cap_seconds:
            print("    (time cap hit, partial)")
            break
        fixed = {}
        for j in range(M):
            for t in range(D):
                fixed[U[j][t]] = 1 if t == dg[j] else 0
        comps = zero_states_clamped(Q, dict(fixed))
        for x in comps:
            assert Q.energy(x) == 0
        is_sol = (dg_point(list(dg)) == Tp)
        collide = not chain_ok(dg)  # degenerate chain (x1==x2 somewhere / hits O)
        status = 'SOL' if is_sol else ('degen' if collide else 'nonsol')
        flag = ''
        if (len(comps) > 0) != is_sol:
            flag = '   *** FAITHFULNESS BREAK ***'
            breaks.append((dg, len(comps), is_sol, collide))
        print(f"    dg={dg} {status:6s} completions={len(comps)}{flag}")
    return breaks


if __name__ == '__main__':
    sys.setrecursionlimit(1000000)
    allbreaks = []
    CASES = [
        (127, 3, 2, 1, 'wallace', 'schoolbook'),
        (127, 3, 2, 1, 'binary', 'schoolbook'),
        (251, 1, 2, 1, 'wallace', 'schoolbook'),
        (127, 3, 3, 1, 'wallace', 'schoolbook'),
        (127, 3, 4, 2, 'wallace', 'schoolbook'),
    ]
    for c in CASES:
        try:
            allbreaks += run(*c)
        except Exception as e:
            import traceback; print(f"  case {c} raised: {e}"); traceback.print_exc()
    print("\nTOTAL FAITHFULNESS BREAKS:", len(allbreaks))
    for b in allbreaks:
        print("  ", b)
