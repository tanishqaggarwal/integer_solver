#!/usr/bin/env python3
"""demo_win.py -- faithfulness of the windowed encoding + window-width scan."""
import time
from ecsmall import curve, find
from ladder import build_win

def make(p, B, m, w, k):
    add, mul = curve(p, B)
    G, order = find(p, B)
    assert order >= (1 << m)
    M = (m + w - 1) // w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(1 << w)] for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))
    T = mul(k, G)
    Tp = add(mul(k % order, G), mul(off % order, G))     # T + offset*G
    return G, order, T, Tp, table, M, add, mul

def chain_ok(table, M, dg, add, mul):
    S = table[0][dg[0]]
    for j in range(1, M):
        Q = table[j][dg[j]]
        if S is None or Q is None or S[0] == Q[0]: return False
        S = add(S, Q)
    return S is not None


def run(p, B, m, w, k, mode='wallace', neq=False, exhaustive=True):
    G, order, T, Tp, table, M, add, mul = make(p, B, m, w, k)
    dgs = lambda kk: [(kk >> (w * j)) % (1 << w) for j in range(M)]
    if not chain_ok(table, M, dgs(k), add, mul):
        k = next(kk for kk in range(1 << m) if chain_ok(table, M, dgs(kk), add, mul))
        G, order, T, Tp, table, M, add, mul = make(p, B, m, w, k)
    t0 = time.time()
    L, U = build_win(p, B, table, Tp, w, mode=mode, neq=neq)
    Q = L.qb; st = Q.stats()
    print(f"[w={w}{' neq' if neq else ''}] p={p} m={m} windows={M}: {st['vars']} vars, {st['couplers']} couplers, "
          f"AND={st['and_vars']}, |J| 2^{st['dynamic_range_bits']}, {time.time()-t0:.2f}s")

    def digits(kk):
        d, x = [], kk
        for j in range(M): d.append(x % (1 << w)); x >>= w
        return d

    def energy_of(kk):
        dg = digits(kk)
        wv0 = {f"_u{j}": dg[j] for j in range(M)}
        inp = {}
        for j in range(M):
            for t in range(1 << w): inp[U[j][t]] = 1 if t == dg[j] else 0
        try: x, _ = Q.witness(inp, wv0)
        except Exception as ex: return None
        return Q.energy(x)

    if energy_of(k) != 0:
        # the ladder for this k passes through a degenerate addition (S_i == +-addend);
        # re-randomise by choosing another k.  Probability at 256 bits is ~2^-247.
        for kk in range(1 << m):
            if mul(kk, G) is not None and energy_of(kk) == 0: k = kk; break
    assert energy_of(k) == 0, "true k is not a zero-energy state"
    if exhaustive:
        zeros = [kk for kk in range(1 << m) if energy_of(kk) == 0]
        degen = sum(1 for kk in range(1 << m) if energy_of(kk) is None)
        T = mul(k, G)
        sol = [kk for kk in range(1 << m) if mul(kk, G) == T]
        assert set(zeros) <= set(sol) and k in zeros, (zeros, sol)
        print(f"    zero-energy k: {zeros}  true solutions: {sol}  "
              f"degenerate-chain k: {degen}/{1<<m}  -- faithful")
    return st

if __name__ == '__main__':
    for w in (1, 2, 3, 4):
        run(97, 3, 4, w, 9)
    print()
    for w in (1, 2, 3, 4):
        run(97, 3, 4, w, 9, neq=True)
