#!/usr/bin/env python3
"""semaev_signs.py -- how many spurious ground states the S_3 sign ambiguity buys,
measured on real Hamiltonians, for BOTH digit conventions.

S_3(x1,x2,x3) = 0  iff  x3 in { x(P+Q), x(P-Q) }.  A chain of M-1 of them
therefore admits 2^{M-1} branch patterns per digit tuple instead of one.  Whether
that costs anything depends entirely on which comb the x-only chain sits on:

  UNSIGNED comb (digits 1..2^w, table (t+1)*2^{wj}G):
      the reference problem has NO sign freedom, so all 2^{M-1} branches beyond
      the intended one are spurious.  Blow-up ~ 2^{M-1}.
  SIGNED comb (digits +-{1,3,..}, table (2t+1)*2^{wj}G):
      the reference problem has exactly one sign per window, i.e. 2^M patterns
      of which 2 give the same point.  The branch freedom and the digit-sign
      freedom are THE SAME freedom, so nothing is lost -- x-only is EXACT.

Both statements are checked here by enumerating every digit tuple and every
branch pattern against the built Hamiltonian.
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ecsmall import curve, find            # noqa: E402
from enc import build_semaev               # noqa: E402


def run(p, B, w, M, signed, mode='wallace', onehot='tree'):
    add, mul = curve(p, B)
    G, order = find(p, B)
    if signed:
        D = 1 << (w - 1)
        pts = [[mul(((2 * t + 1) << (w * j)) % order, G) for t in range(D)]
               for j in range(M)]
    else:
        D = 1 << w
        pts = [[mul(((t + 1) << (w * j)) % order, G) for t in range(D)]
               for j in range(M)]
    table = [[pts[j][t][0] for t in range(D)] for j in range(M)]

    def chain(mags, signs):
        S = pts[0][mags[0]]
        S = (S[0], (-S[1]) % p) if signs[0] else S
        zs = []
        for j in range(1, M):
            P = pts[j][mags[j]]
            P = (P[0], (-P[1]) % p) if signs[j] else P
            if S is None or P is None or S[0] == P[0]:
                return None, None
            S = add(S, P)
            if S is None:
                return None, None
            zs.append(S[0])
        return S, zs

    # target = the point the all-plus reference tuple reaches
    ref = None
    for mags in itertools.product(range(D), repeat=M):
        S, _ = chain(list(mags), [0] * M)
        if S is not None:
            ref = list(mags)
            break
    assert ref is not None
    Sfin, _ = chain(ref, [0] * M)
    xT = Sfin[0]

    L, SEL = build_semaev(p, B, table, xT, w, mode=mode, onehot=onehot)
    Q = L.qb
    st = Q.stats()

    def energy(mags, signs):
        S, zs = chain(mags, signs)
        if S is None:
            return None
        wv = {f"_u{j}": mags[j] for j in range(M)}
        for j in range(1, M - 1):
            wv[f"_z{j}"] = zs[j - 1]
        try:
            x, _ = Q.witness({}, wv)
        except Exception:
            return None
        return Q.energy(x)

    zeros, genuine = [], []
    for mags in itertools.product(range(D), repeat=M):
        ok = False
        for signs in itertools.product((0, 1), repeat=M):
            if energy(list(mags), list(signs)) == 0:
                ok = True
                break
        if not ok:
            continue
        zeros.append(mags)
        # is this digit tuple a solution of the ORIGINAL (y-carrying) problem?
        if signed:
            hit = any(chain(list(mags), list(sg))[0] == Sfin
                      for sg in itertools.product((0, 1), repeat=M))
        else:
            hit = (chain(list(mags), [0] * M)[0] == Sfin)
        if hit:
            genuine.append(mags)
    kind = "SIGNED  " if signed else "UNSIGNED"
    print(f"[{kind} comb, x-only S_3 chain] p={p} w={w} M={M} D={D}: "
          f"{st['vars']} vars")
    print(f"    digit tuples {D ** M}: zero-energy {len(zeros)}, "
          f"genuine {len(genuine)}, SPURIOUS {len(zeros) - len(genuine)}"
          f"   (blow-up {len(zeros) / max(1, len(genuine)):.2f}x, "
          f"2^(M-1) = {2 ** (M - 1)})")
    return len(zeros), len(genuine)


if __name__ == '__main__':
    for (p, B) in ((127, 3), (331, 2)):
        for M in (2, 3, 4):
            try:
                run(p, B, 2, M, signed=False)
                run(p, B, 2, M, signed=True)
            except AssertionError as ex:
                print("  skipped:", ex)
        print()
