#!/usr/bin/env python3
"""ATTACK 1: full ground-state enumeration of a small windowed-comb LADDER.

demo_win2.py only replays ONE canonical witness per scalar kk -- it cannot rule
out a spurious ancilla completion reaching E=0 for a WRONG scalar (or for a
scalar whose canonical completion failed). zero_states() (validated against
true 2^n brute force elsewhere) enumerates the ENTIRE zero-energy set of the
Hamiltonian. We run it on the real ladder2.build_win2 and on ladder.build_win,
then decode every zero-energy state to its digits and check it is a true
discrete-log solution.  Any zero-energy state that is NOT a solution is a
faithfulness break.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.join(HERE, '..', '..', 'squeeze')
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
sys.path.insert(0, os.path.abspath(SQ))
from ecsmall import curve, find
import verify  # for zero_states
from ladder2 import build_win2


def enumerate_ladder(p, B, m, w, mode='wallace', mult='schoolbook', leaf=3):
    add, mul = curve(p, B)
    G, order = find(p, B)
    assert order >= (1 << m), f"order {order} too small for m={m}"
    M = (m + w - 1) // w
    D = 1 << w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(D)] for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))

    # pick a target from a random-ish scalar k whose chain is non-degenerate
    def chain_ok(dg):
        S = table[0][dg[0]]
        for j in range(1, M):
            Q = table[j][dg[j]]
            if S is None or Q is None or S[0] == Q[0]:
                return False
            S = add(S, Q)
        return S is not None
    dgs = lambda kk: [(kk >> (w * j)) % D for j in range((m + w - 1) // w)]
    k = next(kk for kk in range(1 << m) if chain_ok(dgs(kk)))
    Tp = add(mul(k % order, G), mul(off % order, G))

    Q, U = build_win2(p, table, Tp, w, mode=mode, mult=mult, leaf=leaf, red='naf')
    Q.finalize()
    print(f"  p={p} B={B} m={m} w={w} M={M} mode={mode} mult={mult}: vars={Q.n}, k={k}")

    # true solution set: which digit-vectors give sum == Tp
    def digits_to_point(dg):
        S = None
        for j in range(M):
            pt = table[j][dg[j]]
            S = add(S, pt)
        return S
    import itertools
    truth = set()
    for dg in itertools.product(range(D), repeat=M):
        if digits_to_point(list(dg)) == Tp:
            truth.add(tuple(dg))

    t0 = time.time()
    zeros = verify.zero_states(Q)
    dt = time.time() - t0
    # decode each zero-energy state to its digit vector via the one-hot selectors U
    decoded = set()
    bad = []
    for x in zeros:
        assert Q.energy(x) == 0, "zero_states returned a non-zero-energy state!"
        dg = []
        onehot_ok = True
        for j in range(M):
            sel = [t for t in range(D) if x[U[j][t]] == 1]
            if len(sel) != 1:
                onehot_ok = False
                dg.append(('MULTI/NONE', [x[U[j][t]] for t in range(D)]))
            else:
                dg.append(sel[0])
        if not onehot_ok:
            bad.append(('one-hot violated but E=0', dg, x))
            continue
        dgt = tuple(dg)
        decoded.add(dgt)
        if dgt not in truth:
            bad.append(('E=0 decodes to NON-SOLUTION digits', dgt, digits_to_point(list(dgt)), 'want', Tp))
    print(f"    zero-energy states: {len(zeros)}  distinct decoded digit-vecs: {len(decoded)}  "
          f"true solutions: {len(truth)}  ({dt:.1f}s)")
    print(f"    decoded == truth ? {decoded == truth}   spurious(decoded-truth)={sorted(decoded-truth)}   "
          f"missing(truth-decoded)={sorted(truth-decoded)}")
    if bad:
        print("    *** POTENTIAL BREAKS ***")
        for b in bad[:10]:
            print("      ", str(b)[:200])
    return decoded, truth, bad


if __name__ == '__main__':
    sys.setrecursionlimit(1000000)
    allbad = []
    # keep instances tiny enough that zero_states terminates
    CASES = [
        (127, 3, 2, 1, 'wallace', 'schoolbook'),
        (127, 3, 2, 1, 'binary', 'schoolbook'),
        (251, 1, 2, 1, 'wallace', 'schoolbook'),
    ]
    for (p, B, m, w, mode, mult) in CASES:
        try:
            d, t, bad = enumerate_ladder(p, B, m, w, mode=mode, mult=mult)
            allbad += bad
        except Exception as e:
            import traceback
            print(f"  case {(p,B,m,w,mode,mult)} raised: {e}")
            traceback.print_exc()
    print()
    print("TOTAL POTENTIAL BREAKS:", len(allbad))
