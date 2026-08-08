#!/usr/bin/env python3
"""sa_probe.py -- is the annealing landscape of this encoding workable at all?

  P0  sanity: SA on random QUBOs vs brute force  (is the annealer itself correct?)
  P1  ancillas only: clamp the answer digits to the KNOWN solution and ask SA to fill
      in the deterministic ancillas.  Easiest possible task on this Hamiltonian:
      a unique forced completion, zero search over k.
  P2  the real thing: nothing clamped.
"""
import random, time
from ecsmall import curve, find
from ladder import build_win
from sa import SA

P, B = 1021, 3


def p0():
    rnd = random.Random(7); bad = 0
    for trial in range(5):
        n = 12; Q = {}
        for i in range(n):
            Q[(i,)] = rnd.randrange(-9, 10)
            for j in range(i + 1, n):
                if rnd.random() < .4: Q[(i, j)] = rnd.randrange(-9, 10)
        e, _ = SA(Q, n).run(sweeps=400, seed=trial)
        best = min(sum(c * (all((v >> i & 1) for i in m) if m else 1) for m, c in Q.items())
                   for v in range(1 << n))
        bad += (e != best)
    print(f"P0 sanity: SA matched brute force on {5-bad}/5 random QUBOs", flush=True)


def build(m, w, seed=0):
    add, mul = curve(P, B)
    G, order = find(P, B)
    M = (m + w - 1) // w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(1 << w)] for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))
    rnd = random.Random(seed)
    for _ in range(500):
        k = rnd.randrange(1 << m)
        dg = [(k >> (w * j)) % (1 << w) for j in range(M)]
        S = table[0][dg[0]]; ok = True
        for j in range(1, M):
            Qp = table[j][dg[j]]
            if S is None or S[0] == Qp[0]: ok = False; break
            S = add(S, Qp)
        if ok and S is not None: break
    Tp = add(mul(k, G), mul(off % order, G))
    L, U = build_win(P, B, table, Tp, w, mode='wallace')
    return k, L.qb, U, M, dg


if __name__ == '__main__':
    p0()
    for m, w in ((4, 2), (6, 2)):
        k, Q, U, M, dg = build(m, w)
        sa = SA(Q.Q, Q.n)
        clamp = [U[j][t] for j in range(M) for t in range(1 << w)]
        cs = set(clamp)
        x0 = [0] * Q.n
        for j in range(M):
            for t in range(1 << w): x0[U[j][t]] = 1 if t == dg[j] else 0
        print(f"\ninstance m={m} w={w}: {Q.n} vars, "
              f"{sum(1 for mm in Q.Q if len(mm)==2)} couplers, k={k}", flush=True)
        for sweeps in (2000, 20000, 100000):
            t0 = time.time(); hits = 0; bestE = 10 ** 18
            for r in range(6):
                rr = random.Random(r)
                start = [x0[i] if i in cs else rr.randrange(2) for i in range(Q.n)]
                e, _ = sa.run(sweeps=sweeps, seed=r, x0=start, clamp=clamp)
                bestE = min(bestE, e); hits += (e == 0)
            print(f"  P1 answer-clamped, {sweeps:6d} sweeps: E=0 in {hits}/6 runs "
                  f"(best E={bestE})  [{time.time()-t0:.0f}s]", flush=True)
        t0 = time.time(); hits = 0; bestE = 10 ** 18
        for r in range(6):
            e, _ = sa.run(sweeps=20000, seed=100 + r)
            bestE = min(bestE, e); hits += (e == 0)
        print(f"  P2 free search,     20000 sweeps: E=0 in {hits}/6 runs "
              f"(best E={bestE})  [{time.time()-t0:.0f}s]", flush=True)
