#!/usr/bin/env python3
"""subsetsum.py -- the encoding this problem WANTS to be, and what it costs.

E(F_p) is cyclic of prime order n, so as a group it IS Z_n.  Under that
isomorphism the whole instance collapses to a modular subset-sum:

        find c_j in {0,1}  with  sum_j c_j * l_j  ==  l_T   (mod n)

where l_j = dlog_G(Q_j).  That QUBO is tiny -- no field arithmetic at all, one
variable per point plus a quotient word.  We measure it here.

The catch is exactly one thing: you need the l_j.  For the 256 points we HAVE
them (l_i = 2^i, that is what the doubling chain means) -- but the one log we
need, l_T = k, is the answer.  Any set of extra relations with known logs would
drop the encoding from ~9e6 qubits to the number printed below; producing such
relations for a prime-field curve is the open problem, not the encoding.
"""
import random, math
from qubo import QB


def build(logs, target, n, mode='wallace', chunk=16):
    """QUBO for  sum c_j*logs[j] == target (mod n)."""
    Q = QB(mode=mode, chunk=chunk)
    m = len(logs)
    c = [Q.new(f"c{j}", 'input') for j in range(m)]
    for j in range(m):
        Q.trace.append(('word', f"c{j}", [c[j]], (lambda wv, j=j: wv[f"_c{j}"])))
    qlo, qhi = (-target) // n, (sum(logs) - target) // n     # the quotient is signed
    nb = max(1, (qhi - qlo).bit_length())
    q = Q.word("q", nb, lambda wv, qlo=qlo: (sum(logs[j] * wv[f"_c{j}"]
                                                 for j in range(len(logs))) - target) // n - qlo)
    poly = {(c[j],): logs[j] for j in range(m)}
    for t, v in enumerate(q): poly[(v,)] = poly.get((v,), 0) - n * (1 << t)
    Q.assert_zero(poly, -target - n * qlo, "ss")
    Q.finalize()
    return Q, c


def faithful(m=12, bits=20, seed=1):
    rnd = random.Random(seed)
    n = 1048573                                  # prime
    logs = [rnd.randrange(n) for _ in range(m)]
    truth = [rnd.randrange(2) for _ in range(m)]
    target = sum(l * b for l, b in zip(logs, truth)) % n
    Q, c = build(logs, target, n)
    zeros = []
    for v in range(1 << m):
        bits_ = [(v >> j) & 1 for j in range(m)]
        if sum(l * b for l, b in zip(logs, bits_)) % n != target: continue
        x, _ = Q.witness({c[j]: bits_[j] for j in range(m)},
                         {f"_c{j}": bits_[j] for j in range(m)})
        zeros.append(Q.energy(x))
    # and a sample of non-solutions
    bad = 0
    for _ in range(300):
        v = rnd.randrange(1 << m)
        bits_ = [(v >> j) & 1 for j in range(m)]
        if sum(l * b for l, b in zip(logs, bits_)) % n == target: continue
        x, _ = Q.witness({c[j]: bits_[j] for j in range(m)},
                         {f"_c{j}": bits_[j] for j in range(m)})
        bad += (Q.energy(x) == 0)
    print(f"faithfulness (m={m}, n={n}): every true solution has E=0: "
          f"{all(e == 0 for e in zeros)} ({len(zeros)} solutions); "
          f"spurious zero-energy non-solutions in 300 samples: {bad}")


if __name__ == '__main__':
    faithful()
    N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    rnd = random.Random(7)
    print(f"\n{'relations m':>12} {'mode':>8} {'qubits':>9} {'couplers':>10} {'|J|':>7}")
    for m in (256, 512, 1024):
        for mode in ('binary', 'wallace'):
            logs = [rnd.randrange(N) for _ in range(m)]
            Q, _ = build(logs, rnd.randrange(N), N, mode=mode)
            st = Q.stats()
            print(f"{m:12d} {mode:>8} {st['vars']:9,d} {st['couplers']:10,d} "
                  f"{'2^%d' % st['dynamic_range_bits']:>7}")
    print("""
So: with relations in hand the whole instance is a few thousand qubits -- it
fits a real annealer today, with room to spare.  Without them the same decision
problem costs 9.06e6.  The entire 3-orders-of-magnitude gap is the price of one
missing thing: a way to produce elliptic-curve relations with known logarithms
over a PRIME field.  That is the open problem; the encoding is not.""")
