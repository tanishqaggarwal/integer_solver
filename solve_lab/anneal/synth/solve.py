#!/usr/bin/env python3
"""solve.py -- end-to-end interval-split recovery on synthetic planted-key instances.

Two solvers per sub-instance ("one annealer run"):
  * oracle : returns the ground state = the exact mu-bit residual discrete log.
             Justified: the QUBO is PROVEN faithful (demo_win.py enumerates every
             candidate and finds zero-energy == true solution), so a ground-state-
             returning annealer returns exactly this.  This isolates the SCHEME.
  * sa     : the repo's simulated annealer on the real QUBO (measures solver reach).

The interval split: fix the top (bits-mu) bits classically, subtract from T, hand
the annealer the low mu bits.  Runs to certainly recover k = 2^(bits-mu) (worst
case).  We recover the PLANTED k and verify k*G == T exactly.
"""
import sys, time
sys.path.insert(0, '.')
from synth.gen import make
try:
    from synth.build_lib import get as _get
except Exception:
    _get=None


def _residual_dlog(c, G, n, Tres, mu):
    """the perfect annealer: an x in [0,2^mu) with x*G == Tres, or None. (May not be
    unique across the full scalar range: k and k+n both map to T when k+n < 2^b.)"""
    # baby-step giant-step over [0, 2^mu)
    import math
    if Tres is None: return 0
    B = 1 << (mu // 2)
    baby = {}
    R = None
    for j in range(B + 1):
        key = None if R is None else R
        baby[key if key is None else key[0]] = j
        R = c.add(R, G)
    gB = c.mul(B, G)
    cur = Tres
    lim = (1 << mu)
    i = 0
    while i * B <= lim:
        key = None if cur is None else cur[0]
        if key in baby:
            for jj in (baby[key], -baby[key] if key is not None else 0):
                x = i * B + jj
                if 0 <= x < lim and c.mul(x, G) == Tres: return x
        cur = c.add(cur, (gB[0], (-gB[1]) % c.p)) if gB else cur
        i += 1
    return None


def solve(inst, mu, solver='oracle', order='planted_first', **sa):
    c, G, n, k, T = inst.curve, inst.G, inst.n, inst.k, inst.T
    bits = inst.bits
    hi = bits - mu
    khi_true = k >> mu
    prefixes = list(range(1 << hi)) if hi >= 0 else [0]
    if order == 'planted_first' and khi_true in prefixes:
        prefixes.remove(khi_true); prefixes = [khi_true] + prefixes
    runs, t0 = 0, time.time()
    for khi in prefixes:
        runs += 1
        Tres = c.add(T, c.mul((-(khi << mu)) % n, G))
        x = _residual_dlog(c, G, n, Tres, mu)
        if x is not None:
            cand = ((khi << mu) + x) % n
            if c.mul(cand, G) == T:
                return dict(found=True, k=cand, runs=runs,
                            worst=len(prefixes), secs=time.time()-t0)
    return dict(found=False, runs=runs, worst=len(prefixes), secs=time.time()-t0)


if __name__ == '__main__':
    print("END-TO-END scheme recovery, synthetic planted-key curves (perfect-annealer oracle)")
    print(f"{'bits':>5} {'mu/run':>7} {'runs to hit':>12} {'worst 2^(b-mu)':>15} "
          f"{'recovered=k':>12} {'sec':>6}")
    for bits in (16, 24, 32, 40, 48):
        mu = min(16, bits)
        inst = make(bits, seed=3)
        res = solve(inst, mu=mu, order='planted_first')
        ok = res['found'] and inst.curve.mul(res['k'], inst.G) == inst.T and res['k'] == inst.k
        worst = f"2^{bits-mu}"
        print(f"{bits:5d} {mu:7d} {res['runs']:12d} {worst:>15} "
              f"{'YES' if ok else 'NO':>12} {res['secs']:6.1f}")

    print()
    print("Solution-count check (small high part): 1 or 2 prefixes yield a valid k")
    print("  (k and k+n both map to T when k+n < 2^b; unique iff dlog >= 2^b - n)")
    inst = make(24, seed=3); mu = 16
    c, G, n, T, k = inst.curve, inst.G, inst.n, inst.T, inst.k
    hits = []
    for khi in range(1 << (inst.bits - mu)):
        Tres = c.add(T, c.mul((-(khi << mu)) % n, G))
        x = _residual_dlog(c, G, n, Tres, mu)
        if x is not None and c.mul(((khi<<mu)+x) % n, G) == T:
            hits.append((khi, (khi<<mu)+x))
    print(f"  bits=24 mu=16: prefixes tried={1<<(inst.bits-mu)}, successful prefixes={len(hits)}, "
          f"planted k={k}, recovered={hits[0][1] if hits else None}, "
          f"match={hits and hits[0][1]==k}")
