#!/usr/bin/env python3
"""unit_probe.py -- how big a piece of this Hamiltonian can an annealer actually settle?

The atom of the encoding is one modular multiplication  a*b == c (mod p).
Two questions, both measured with plain simulated annealing:
  Q1  with a and b CLAMPED to a known pair, can the annealer fill in the forced
      ancillas (partial products, adder bits, quotient)?  If not, no amount of
      qubits helps -- the arithmetic layer itself is not annealable.
  Q2  with nothing clamped, can it invert the multiplication?
"""
import random, time, json
from ladder import Ladder
from sa import SA


def modmul_qubo(s, mode='wallace', seed=3):
    rnd = random.Random(seed)
    p = (1 << (s - 1)) + 2 * rnd.randrange(1 << (s - 3)) + 1
    while True:
        ok = all(p % q for q in range(3, int(p ** .5) + 1, 2))
        if ok: break
        p += 2
    a, b = rnd.randrange(2, p), rnd.randrange(2, p)
    c = a * b % p
    L = Ladder(p, mode=mode); Q = L.qb
    A = Q.word("a", s, lambda wv: a)
    Bw = Q.word("b", s, lambda wv: b)
    L.mul_eq("mm", A, Bw, "a", "b", [], c)
    Q.finalize()
    return p, a, b, c, Q, A, Bw


if __name__ == '__main__':
    print(f"{'s':>3} {'p':>8} {'vars':>7} {'coupl':>8} | {'Q1 clamped':>22} | {'Q2 free':>18}")
    res = []
    for s in (6, 8, 10, 12):
        p, a, b, c, Q, A, Bw = modmul_qubo(s)
        x, _ = Q.witness({}, {})
        assert Q.energy(x) == 0
        sa = SA(Q.Q, Q.n)
        clamp = A + Bw
        x0 = [0] * Q.n
        for t, v in enumerate(A): x0[v] = (a >> t) & 1
        for t, v in enumerate(Bw): x0[v] = (b >> t) & 1
        h1 = b1 = None
        t0 = time.time(); hits = 0; best = 10**18
        for r in range(8):
            rr = random.Random(r)
            st = [x0[i] if i in set(clamp) else rr.randrange(2) for i in range(Q.n)]
            e, _ = sa.run(sweeps=8000, seed=r, x0=st, clamp=clamp)
            best = min(best, e); hits += (e == 0)
        h1, b1, t1 = hits, best, time.time() - t0
        t0 = time.time(); hits = 0; best = 10**18
        for r in range(8):
            e, _ = sa.run(sweeps=8000, seed=50 + r)
            best = min(best, e); hits += (e == 0)
        print(f"{s:3d} {p:8d} {Q.n:7d} {sum(1 for m in Q.Q if len(m)==2):8d} | "
              f"{h1}/8 hit, best E={b1:<6} {t1:5.0f}s | {hits}/8 hit, best E={best}")
        res.append(dict(s=s, p=p, vars=Q.n, clamped_hits=h1, clamped_best=b1,
                        free_hits=hits, free_best=best))
    json.dump(res, open('unit_probe.json', 'w'), indent=1)
