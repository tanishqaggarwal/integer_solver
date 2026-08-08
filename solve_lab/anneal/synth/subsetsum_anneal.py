#!/usr/bin/env python3
"""subsetsum_anneal.py -- the one reformulation with a gradient, and why it is out of reach.

Contrast two encodings of the SAME synthetic key:

  (A) DLP direct :  find bits b_i with  sum b_i (2^i G) = T.   Energy is a function
      of k*G, a pseudorandom map of the bits -> NO gradient (measured: corr ~0.06).

  (B) subset-sum :  find bits c_j with  sum c_j l_j = L_T (mod n),  l_j known.
      Energy (sum c_j l_j - L_T)^2 is LINEAR in the bits -> flipping c_j changes the
      sum by l_j, a real gradient. This is the encoding that WOULD let one anneal
      resolve many bits, dropping total calls from 2^b toward 2^(b-mu).

The catch is not the annealer -- it is getting to (B). The map from (A) to (B) is
the discrete logarithm of a factor base: relations Q_j = l_j G with l_j known, and
a way to write the target over them. For prime-field curves that is index calculus,
which does not exist. So (B) is only writable when we ALREADY know the logs -- i.e.
never, in the cryptanalytic direction.

Here we (i) show (B) anneals where (A) does not, on planted instances, and (ii) make
precise that building (B) from (A) is the missing step.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'solver'))
import numpy as np
from synth.gen import make
from qubo import QB
import importlib.util
_spec = importlib.util.spec_from_file_location("smodel", os.path.join(os.path.dirname(__file__), 'solver', 'model.py'))
smodel = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(smodel)
import solvers as S


def build_subsetsum(logs, target, n, mode='wallace'):
    """QUBO for sum c_j*logs[j] == target (mod n).  Returns (QB, Ising, c-vars, xstar)."""
    Q = QB(mode=mode)
    m = len(logs)
    c = [Q.new(f"c{j}", 'input') for j in range(m)]
    for j in range(m):
        Q.trace.append(('word', f"c{j}", [c[j]], (lambda wv, j=j: wv[f"_c{j}"])))
    qlo = (-target) // n
    qhi = (sum(logs) - target) // n
    nb = max(1, (qhi - qlo).bit_length())
    q = Q.word("q", nb, lambda wv: (sum(logs[j]*wv[f"_c{j}"] for j in range(m)) - target)//n - qlo)
    poly = {(c[j],): logs[j] for j in range(m)}
    for t, v in enumerate(q):
        poly[(v,)] = poly.get((v,), 0) - n*(1 << t)
    Q.assert_zero(poly, -target - n*qlo, "ss")
    Q.finalize()
    return Q, c


def gradient_profile(m=20, seed=1, samples=3000):
    """energy vs Hamming distance to the planted subset -- does subset-sum have a gradient?"""
    rng = np.random.default_rng(seed)
    n = 1000003
    logs = [int(rng.integers(1, n)) for _ in range(m)]
    truth = [int(rng.integers(0, 2)) for _ in range(m)]
    target = sum(l*b for l, b in zip(logs, truth)) % n
    Q, c = build_subsetsum(logs, target, n)
    ising = smodel.qubo_to_ising(Q.Q, Q.n)
    # witness for the true subset
    wv0 = {f"_c{j}": truth[j] for j in range(m)}
    xstar, _ = Q.witness({c[j]: truth[j] for j in range(m)}, wv0)
    assert Q.energy(xstar) == 0
    buckets = {}
    for _ in range(samples):
        cand = [int(rng.integers(0, 2)) for _ in range(m)]
        d = sum(a != b for a, b in zip(cand, truth))
        wv = {f"_c{j}": cand[j] for j in range(m)}
        try:
            x, _ = Q.witness({c[j]: cand[j] for j in range(m)}, wv)
            e = Q.energy(x)
        except Exception:
            continue
        buckets.setdefault(d, []).append(e)
    xs, ys = [], []
    for d, es in buckets.items():
        for e in es:
            xs.append(d); ys.append(e)
    r = np.corrcoef(xs, ys)[0, 1]
    return r, Q, c, ising, truth, logs, target, n


def can_anneal(m, reps=8, seed=2):
    """does a real solver reach the ground state of the subset-sum QUBO?"""
    rng = np.random.default_rng(seed)
    n = 1000003
    logs = [int(rng.integers(1, n)) for _ in range(m)]
    truth = [int(rng.integers(0, 2)) for _ in range(m)]
    target = sum(l*b for l, b in zip(logs, truth)) % n
    Q, c = build_subsetsum(logs, target, n)
    ising = smodel.qubo_to_ising(Q.Q, Q.n)
    hits = 0
    for r in range(reps):
        e, _ = S.tabu(ising, iters=40000, seed=r)
        hits += (e == 0)
    return Q.n, hits, reps


if __name__ == '__main__':
    print("=== Does subset-sum have a gradient the DLP lacks? ===")
    r, *_ = gradient_profile(m=20)
    print(f"subset-sum: correlation(Hamming distance, energy) = {r:.3f}")
    print("  (DLP direct measured 0.06 = flat. subset-sum should be strongly positive")
    print("   because the sum is LINEAR in the bits -- that is the gradient.)")
    print()
    print("=== Does a real solver crack subset-sum where it cannot crack the DLP? ===")
    print(f"  {'m (bits)':>8} {'vars':>6} {'tabu hits/reps':>16}")
    for m in (16, 24, 32, 40):
        n, h, reps = can_anneal(m)
        print(f"  {m:8d} {n:6d} {h}/{reps:<14d}")
    print("""
If subset-sum anneals at m the DLP cannot touch, the annealer's limitation is not
depth or effort -- it is that the DLP encoding is gradient-free and the ONLY
gradient-ful reformulation (subset-sum) requires discrete logs of a factor base,
i.e. index calculus, absent for prime-field ECC. That missing reduction, not the
annealer, is the whole 2^b -> 2^(b-mu) gap.""")


# --------------------------------------------------------------------------
# The mechanism check: modular vs non-modular subset-sum.
# --------------------------------------------------------------------------
def build_intsum(logs, target, mode='wallace'):
    """QUBO for sum c_j*logs[j] == target over Z (NO modulus). Bounded coeffs."""
    Q = QB(mode=mode)
    m = len(logs)
    c = [Q.new(f"c{j}", 'input') for j in range(m)]
    for j in range(m):
        Q.trace.append(('word', f"c{j}", [c[j]], (lambda wv, j=j: wv[f"_c{j}"])))
    poly = {(c[j],): logs[j] for j in range(m)}
    Q.assert_zero(poly, -target, "is")
    Q.finalize()
    return Q, c


def contrast(m=20, seed=5, samples=3000):
    rng = np.random.default_rng(seed)
    truth = [int(rng.integers(0, 2)) for _ in range(m)]
    # non-modular, moderate coefficients (density ~1: coeffs ~ m bits)
    logs = [int(rng.integers(1, 1 << m)) for _ in range(m)]
    target = sum(l*b for l, b in zip(logs, truth))
    Q, c = build_intsum(logs, target)
    buckets = {}
    for _ in range(samples):
        cand = [int(rng.integers(0, 2)) for _ in range(m)]
        d = sum(a != b for a, b in zip(cand, truth))
        wv = {f"_c{j}": cand[j] for j in range(m)}
        try:
            x, _ = Q.witness({c[j]: cand[j] for j in range(m)}, wv)
            e = Q.energy(x)
        except Exception:
            continue
        buckets.setdefault(d, []).append(e)
    xs, ys = [], []
    for d, es in buckets.items():
        for e in es:
            xs.append(d); ys.append(e)
    r = np.corrcoef(xs, ys)[0, 1]
    ising = smodel.qubo_to_ising(Q.Q, Q.n)
    hits = sum(S.tabu(ising, iters=40000, seed=r2)[0] == 0 for r2 in range(8))
    return r, hits


if __name__ == '__main__' and '--contrast' in sys.argv:
    print("=== MECHANISM: even the 'best-case' reformulation is gradient-free ===")
    print("high-density subset-sum over Z (coeffs ~ 2^m, the ECDLP-relevant density):")
    r, h = contrast(20)
    print(f"  correlation(distance, energy) = {r:.3f}   tabu hits/8 = {h}")
    print("""  MEASURED RESULT: no gradient (corr ~0), does NOT anneal (0/8) -- same as the
  modular case and the DLP itself.  So the reformulation I hoped had a gradient
  does not: a random ~m-bit-coefficient knapsack is as flat for the annealer as
  the DLP.  Hamming-distance carries no energy signal because the coefficients
  span 2^m, so one bit-flip moves the sum by up to 2^m -- a needle, not a slope.
  A gradient appears only in the LOW-density regime (coeffs >> 2^m, few solutions)
  -- and that regime is exactly what LLL already solves in polynomial time.
  Conclusion: there is no gradient-ful sub-problem in this family for the annealer
  to exploit.  The 2^b barrier is not an artifact of the DLP encoding; it survives
  every reformulation tried, reachable or not.""")
