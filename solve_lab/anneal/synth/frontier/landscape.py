#!/usr/bin/env python3
"""landscape.py -- diagnose.py's E-vs-Hamming-distance profile, ported to the
MODMUL.  Distance is measured in the FREE OPERAND bits (b): operand a is pinned
to its planted value, operand b is set to a candidate, and every ancilla is
FORCED to the value the witness gives for (a, candidate_b) -- exactly as
diagnose.py forces the comb's ancillas.  E then reflects only the arithmetic
verifier's signal: does a wrong b sit LOWER (a gradient a solver can follow) or
just as high as a random b (a flat plateau -- a needle)?

We do NOT modify model.py: we rebuild the same a*b==c(mod p) QUBO here with b's
witness function reading a candidate from wv0, so Q.witness can fill the ancillas
for any candidate b.
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'synth', 'solver'))
import numpy as np
from ladder import Ladder
import fbuild as FB


def build_probe(s, mode='wallace', seed=3):
    """a*b==c (mod p) with b's value read from wv0['_bcand'] so any candidate b
    can be witnessed.  Planted (a,b,c) identical to model.build_modmul."""
    p, a, b, cc = FB.pick_prime(s, seed=seed)
    L = Ladder(p, mode=mode)
    Q = L.qb
    A = Q.word("a", s, lambda wv, a=a: a)
    Bw = Q.word("b", s, lambda wv, b=b: wv.get("_bcand", b))
    L.mul_eq("mm", A, Bw, "a", "b", [], cc)
    Q.finalize()
    # sanity: planted candidate is E=0
    x, _ = Q.witness({}, {"_bcand": b})
    assert Q.energy(x) == 0
    return dict(Q=Q, p=p, a=a, b=b, c=cc, s=s, Bbits=Bw)


def profile(s, samples=4000, seed=1, mode='wallace'):
    pr = build_probe(s, mode=mode)
    Q, p, a, b = pr['Q'], pr['p'], pr['a'], pr['b']
    rng = np.random.default_rng(seed)
    buckets = {}
    # enumerate all b in [0,p) if small, else sample
    cands = range(p) if p <= samples else (int(rng.integers(0, p)) for _ in range(samples))
    for cand in cands:
        d = bin((cand ^ b) & ((1 << s) - 1)).count('1')
        try:
            x, _ = Q.witness({}, {"_bcand": cand})
            e = Q.energy(x)
        except Exception:
            continue
        buckets.setdefault(d, []).append(e)
    rows = []
    xs, ys = [], []
    for d in sorted(buckets):
        es = buckets[d]
        rows.append((d, len(es), min(es), float(np.mean(es)), max(es)))
        for e in es:
            xs.append(d); ys.append(e)
    xs, ys = np.array(xs), np.array(ys)
    r = float(np.corrcoef(xs, ys)[0, 1]) if len(set(xs.tolist())) > 1 else float('nan')
    # gradient test: is min-E at distance 1 meaningfully below min-E at large d?
    return dict(s=s, p=p, rows=rows, corr=r, n_valid=len(xs))


def report(s, mode='wallace', samples=4000):
    pf = profile(s, samples=samples, mode=mode)
    print(f"\nMODMUL s={s} (p={pf['p']}, mode={mode}): E vs distance-to-planted-b "
          f"[{pf['n_valid']} valid candidates]")
    print(f"  (operand a pinned; ancillas FORCED to witness(a,b_cand))")
    print(f"  {'dist':>4} {'n':>5} {'minE':>7} {'meanE':>9} {'maxE':>7}")
    for d, n, mn, mean, mx in pf['rows']:
        print(f"  {d:4d} {n:5d} {mn:7.0f} {mean:9.1f} {mx:7.0f}")
    print(f"  correlation(distance, energy) = {pf['corr']:.3f}")
    return pf


if __name__ == '__main__':
    import json
    out = {}
    for s in [5, 6, 7, 8, 9, 10]:
        pf = report(s)
        out[s] = dict(p=pf['p'], corr=pf['corr'], n_valid=pf['n_valid'],
                      rows=pf['rows'])
    json.dump(out, open(os.path.join(HERE, 'landscape.json'), 'w'), indent=1)
    print("\nLANDSCAPE DONE")
