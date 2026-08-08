#!/usr/bin/env python3
"""restart_scaling.py -- the run-count crux.

If the landscape has no gradient (diagnose.py), a real annealer finds the needle
only by chance, so the per-anneal success probability should fall like 2^-mu and
the ANNEALS-per-sub-instance should rise like 2^mu.  If so, then

    total anneals = (outer runs 2^(b-mu)) x (anneals per sub-instance ~2^mu) ~ 2^b

is INVARIANT in mu: shrinking the encoding to raise mu cuts outer runs but makes
each inner solve exponentially harder, and the annealer yields no speedup over
brute force.  Measured here by counting SA restarts to first E=0.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath('.'))
from synth.gen import make
import model as M, solvers as S
import numpy as np


def hits_per_restart(bits, mu, w=2, sweeps=3000, restarts=400, seed=0):
    inst = make(bits, seed=3)
    md = M.build_comb(inst, mu, w=w, mode='wallace')
    ising = md['ising']
    hits = 0
    t0 = time.time()
    for r in range(restarts):
        e, _ = S.sa(ising, sweeps=sweeps, seed=seed * 10007 + r)
        if e == 0:
            hits += 1
    p = hits / restarts
    return md['Q'].n, p, time.time() - t0


if __name__ == '__main__':
    print("SA hit-rate per restart vs mu (needle density). sweeps fixed per restart.")
    print(f"{'bits':>4} {'mu':>3} {'n':>6} {'restarts':>8} {'hits':>5} "
          f"{'p(hit)':>9} {'1/p':>8} {'2^mu':>8}")
    for mu in (2, 4, 6, 8):
        bits = max(12, mu + 4)
        n, p, dt = hits_per_restart(bits, mu, restarts=300, sweeps=3000, seed=mu)
        pr = f"{p:.4f}" if p else "<.0033"
        inv = f"{1/p:.0f}" if p else ">300"
        print(f"{bits:4d} {mu:3d} {n:6d} {300:8d} {int(p*300):5d} {pr:>9} {inv:>8} {1<<mu:8d}")
    print("\nIf 1/p tracks 2^mu, the annealer is random sampling: no speedup over")
    print("brute force, and total anneals ~ 2^bits regardless of how mu is chosen.")
