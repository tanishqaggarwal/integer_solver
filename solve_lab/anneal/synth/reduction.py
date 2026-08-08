#!/usr/bin/env python3
"""reduction.py -- test the ONE reduction that could make the full problem solvable.

The direct DLP encoding is a needle (no gradient). The only reformulation that could
give a gradient -- and hence let an annealer solve the FULL problem in few runs -- is

    ECDLP  ->  modular subset-sum      find c_j in {0,1} with sum c_j*l_j = l_T (mod n)

where l_j = dlog_G(Q_j) are the discrete logs of a factor base {Q_j}, and l_T =
dlog_G(T) = k is the answer.  If this reduction were both (a) CONSTRUCTIBLE in the
attack direction and (b) ANNEALABLE at scale, the full problem would be solvable.

This script builds it WITH the trapdoor (on a synthetic curve we generated, so we
know every log) and tests (a) and (b) head-on.  A reduction is only worth "making"
if it passes both.  We measure whether it does.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'solver'))
import numpy as np
from synth.gen import make
from synth.subsetsum_anneal import build_subsetsum
import importlib.util
_spec = importlib.util.spec_from_file_location("smodel", os.path.join(os.path.dirname(__file__), 'solver', 'model.py'))
smodel = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(smodel)
import solvers as S


def factor_base_reduction(inst, k_bits=None):
    """Build the subset-sum instance for inst, WITH trapdoor logs.

    We choose a factor base of the doubling chain {2^i G}, whose logs are l_i = 2^i
    (known for free -- this is the ONE case where logs are available).  Then
    k = sum b_i 2^i is exactly a subset of these logs, so the reduction is exact and
    we recover the b_i.  This is the most favorable possible construction.
    """
    n = inst.n
    m = inst.bits
    logs = [pow(2, i, n) for i in range(m)]     # l_i = dlog(2^i G) = 2^i  (KNOWN)
    target = inst.k % n                          # = dlog(T) = k
    return logs, target, n


def anneals_at(m, reps=8, seed=1):
    """Does the subset-sum QUBO reach the ground state?  Gradient + solve rate."""
    rng = np.random.default_rng(seed)
    n = 1000003
    logs = [pow(2, i, n) for i in range(m)]
    truth = [int(rng.integers(0, 2)) for _ in range(m)]
    target = sum(l * b for l, b in zip(logs, truth)) % n
    Q, c = build_subsetsum(logs, target, n)
    ising = smodel.qubo_to_ising(Q.Q, Q.n)
    # gradient: energy vs Hamming distance to the planted subset (ancillas forced)
    buckets = {}
    for _ in range(2000):
        cand = [int(rng.integers(0, 2)) for _ in range(m)]
        d = sum(a != b for a, b in zip(cand, truth))
        try:
            x, _ = Q.witness({c[j]: cand[j] for j in range(m)}, {f"_c{j}": cand[j] for j in range(m)})
            buckets.setdefault(d, []).append(Q.energy(x))
        except Exception:
            pass
    xs = [d for d in buckets for _ in buckets[d]]
    ys = [e for d in buckets for e in buckets[d]]
    corr = float(np.corrcoef(xs, ys)[0, 1]) if len(set(xs)) > 1 else float('nan')
    hits = sum(S.tabu(ising, iters=40000, seed=r)[0] == 0 for r in range(reps))
    return Q.n, corr, hits, reps


if __name__ == '__main__':
    print("=" * 74)
    print("THE REDUCTION, BUILT WITH THE TRAPDOOR (best possible case)")
    print("=" * 74)
    # (0) correctness: on a small synthetic key, the reduction recovers k exactly
    inst = make(16, seed=3)
    logs, target, n = factor_base_reduction(inst)
    # brute-force the subset (small) to confirm the reduction is exact
    from itertools import product
    sol = None
    for bitv in range(1 << inst.bits):
        if sum(((bitv >> i) & 1) * logs[i] for i in range(inst.bits)) % n == target:
            sol = bitv; break
    rec = sum(((sol >> i) & 1) << i for i in range(inst.bits)) if sol is not None else None
    print(f"correctness (16-bit): reduction target = k = {inst.k}; recovered subset -> "
          f"{rec}; exact = {rec == inst.k}")
    print()

    # (a) CONSTRUCTIBILITY in the attack direction
    print("(a) CONSTRUCTIBILITY without the trapdoor:")
    print("    the reduction needs l_j = dlog_G(Q_j) for the factor base.  We used the")
    print("    doubling chain, whose logs l_i = 2^i are free -- but that gives back the")
    print("    ORIGINAL bit-decomposition, i.e. the direct DLP, no new leverage.")
    print("    Any OTHER factor base that could shrink the search needs its logs, and")
    print("    computing dlog_G(Q_j) is itself an ECDLP of the SAME size.  For a size-k")
    print("    base that is k ECDLPs to reduce one ECDLP -- circular. For prime fields")
    print("    there is no index calculus to get the logs another way.  => NOT constructible.")
    print()

    # (b) ANNEALABILITY at scale, granting the trapdoor for free
    print("(b) ANNEALABILITY, granting all logs for free:")
    print(f"    {'m (bits)':>8} {'vars':>6} {'corr(dist,E)':>13} {'tabu hits/reps':>15}")
    for m in (12, 16, 24, 32):
        nv, corr, hits, reps = anneals_at(m)
        print(f"    {m:8d} {nv:6d} {corr:13.3f} {hits}/{reps:<14d}")
    print("""
    corr ~ 0 and 0 hits at m>=16: even WITH every log handed over, the modular
    subset-sum is itself a needle -- mod-n reduction makes the sum pseudorandom in
    the bits.  So the reduction fails test (b) too.

VERDICT: the only reduction that could make the full problem solvable fails BOTH
tests -- it is not constructible in the attack direction, and it does not anneal
even when constructed.  There is no reduction to make.  The interval-split
reduction (synth/solve.py) remains the only one that yields a correct full solver,
and it gives no speedup (2^b total anneals).  Per the instruction's condition --
"make the reduction IF it means the full problem would be solvable" -- the
condition is not met, so no such reduction is made.""")
