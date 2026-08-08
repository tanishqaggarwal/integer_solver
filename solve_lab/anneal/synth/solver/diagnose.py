#!/usr/bin/env python3
"""diagnose.py -- is the hardness energy-barriers (annealer helps) or no-gradient
(nothing helps)?  Measure E vs Hamming distance to the unique ground state on a
solvable instance, over the ANSWER digits (ancillas set to their forced values)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath('.'))
from synth.gen import make
import model as M
import numpy as np


def profile(bits=12, mu=8, w=2, samples=3000, seed=1):
    inst = make(bits, seed=3)
    md = M.build_comb(inst, mu, w=w, mode='wallace')
    Q, U, Mw, ww = md['Q'], md['U'], md['M'], md['w']
    rng = np.random.default_rng(seed)
    true_dig = md['digits']
    buckets = {}
    for _ in range(samples):
        dg = [int(rng.integers(0, 1 << ww)) for _ in range(Mw)]
        d = sum((a != b) for a, b in zip(dg, true_dig))
        wv0 = {f"_u{j}": dg[j] for j in range(Mw)}
        inp = {}
        for j in range(Mw):
            for t in range(1 << ww):
                inp[U[j][t]] = 1 if t == dg[j] else 0
        try:
            x, _ = Q.witness(inp, wv0)
            e = Q.energy(x)
        except Exception:
            continue
        buckets.setdefault(d, []).append(e)
    print(f"instance bits={bits} mu={mu} w={w}: E vs digit-distance-to-solution")
    print(f"  (ancillas FORCED to each candidate; E then measures the verifier's signal)")
    print(f"  {'dist':>4} {'n':>5} {'minE':>7} {'meanE':>9} {'maxE':>7}")
    for d in sorted(buckets):
        es = buckets[d]
        print(f"  {d:4d} {len(es):5d} {min(es):7.0f} {np.mean(es):9.1f} {max(es):7.0f}")
    xs, ys = [], []
    for d, es in buckets.items():
        for e in es:
            xs.append(d); ys.append(e)
    xs, ys = np.array(xs), np.array(ys)
    r = np.corrcoef(xs, ys)[0, 1] if len(set(xs.tolist())) > 1 else float('nan')
    print(f"  correlation(distance, energy) = {r:.3f}")
    print(f"  -> near 0 means NO gradient: E at distance 1 is no lower than at max")
    print(f"     distance; the annealer sees a flat plateau + one hole (a needle).")


if __name__ == '__main__':
    profile(12, 8, 2)
    print()
    profile(16, 10, 2)
