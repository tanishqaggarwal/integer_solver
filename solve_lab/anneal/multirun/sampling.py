#!/usr/bin/env python3
"""sampling.py -- family S1: the annealer as a SAMPLER inside a van Oorschot-Wiener
distinguished-point collision search.

The scheme.  Fix a t-bit predicate pi on curve points ("distinguished point": the
low t bits of the x-coordinate equal a chosen constant c).  Build a QUBO whose
free inputs are two scalars a in [0,2^mu) and b in [0,2^nu), which computes
P = a*G + b*T (two combs joined by one extra EC addition) and whose only
terminal penalty is pi(P), i.e. a t-bit congruence instead of the full 256-bit
"P == T".  Every zero-energy state is a pair (a,b) whose point is distinguished.
Draw N such pairs; a repeated point between (a,b) and (a',b') gives

        k = (a - a') * (b' - b)^{-1}  mod n .

Accounting (all of it elementary, none of it charitable to the classical side):

  n = 2^256 group elements; the distinguished set has  n / 2^t = 2^(256-t) points.
  Birthday over that set:            N  =  2^((256-t)/2)  =  2^(128 - t/2)  samples.
  The QUBO's solution manifold has   2^(mu+nu-t)  elements, and you cannot draw
  more distinct samples than exist:  mu + nu  >=  128 + t/2 .

Writing MU = mu + nu (the total free scalar width the QPU must carry), the best
choice is the largest legal t, t = 2(MU - 128), which gives

        N = 2^(128 - t/2) = 2^(256 - MU)          runs
        V = V_comb(MU) + one extra EC addition    qubits per run

-- exactly the D1 interval-split curve, shifted UP by one EC addition.  The
square-root the birthday buys is exactly cancelled by having to carry two
scalars instead of one.  Sampling buys nothing.

Writes multirun/sampling.json.
"""
import json, math, os, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from tradeoff import V_comb, V_min, D, HW      # measured atoms, real p


def V_sample(mode, MU, w):
    """two combs of total width MU, joined by one more EC addition.

    The terminal check is a t-bit congruence rather than a 512-bit one, which is
    strictly cheaper than `final`; we charge the full `final` anyway (a bound in
    the scheme's favour is not needed -- it is 0.01% of the total).
    """
    v = V_comb(mode, MU, w)
    return None if v is None else v + D[mode]['add']['vars']


def V_sample_min(mode, MU):
    best = None
    for w in range(1, 13):
        v = V_sample(mode, MU, w)
        if v is not None and (best is None or v < best[0]):
            best = (v, w)
    return best


if __name__ == '__main__':
    out = {'rows': {}}
    for mode in ('binary', 'wallace'):
        print("=" * 104)
        print(f"S1  ANNEALER-AS-SAMPLER INSIDE vOW COLLISION SEARCH   (mode={mode}, real p)")
        print("=" * 104)
        print(f"{'MU=mu+nu':>9} {'t (DP bits)':>12} {'runs N':>10} {'V_S1':>14} "
              f"{'V_D1 (same runs)':>18} {'S1 penalty':>12} {'fits 4400?':>11}")
        rows = []
        for MU in (128, 129, 136, 144, 160, 176, 192, 208, 224, 240, 254, 256):
            t = 2 * (MU - 128)
            if t > 256:
                continue
            N = 256 - MU
            vs, ws = V_sample_min(mode, MU)
            vd, wd = V_min(mode, MU)
            print(f"{MU:9d} {t:12d} {'2^%d' % N:>10} {vs:14,d} {vd:18,d} "
                  f"{'+%.1f%%' % (100*(vs-vd)/vd):>12} {'YES' if vs <= HW else 'no':>11}")
            rows.append(dict(MU=MU, t=t, log2_runs=N, V_S1=vs, V_D1=vd,
                             fits_4400=(vs <= HW)))
        out['rows'][mode] = rows
        print()

    print("=" * 104)
    print("WHAT SAMPLING WOULD HAVE TO BE ABLE TO DO, PER SAMPLE")
    print("=" * 104)
    r = D['binary']
    MU = 192
    v, w = V_sample_min('binary', MU)
    M = -(-MU // w)
    print(f"  Take MU = {MU} (t = {2*(MU-128)} distinguished bits, N = 2^{256-MU} samples).")
    print(f"  Each sample is one anneal of a {v:,d}-qubit QUBO ({M} comb windows, "
          f"{M} EC additions),")
    print(f"  and the anneal must land in the ground-state manifold "
          f"(2^{MU - 2*(MU-128)} states out of 2^{MU} inputs), uniformly,")
    print(f"  with no repeats beyond the birthday budget.  Section 6 of ENCODING.md "
          f"measured that")
    print(f"  simulated annealing cannot reach E=0 even on ONE modular multiplication "
          f"({r['mul']['vars']:,d} qubits)")
    print(f"  with both operands clamped, so the sampling assumption is not merely "
          f"unproven, it is")
    print(f"  contradicted at 1/{v // r['mul']['vars']}th of the required scale.")

    print()
    print("  Also: a QUBO ground state is one assignment. To get a COLLISION inside a")
    print("  single run instead, the QUBO must contain both members of the pair, i.e.")
    print("  4 scalars (a,b,a',b') -- at 256 bits each that is ~4x the monolithic QUBO:")
    for mode in ('binary', 'wallace'):
        v256, _ = V_min(mode, 256)
        print(f"    {mode:>7}: {4*v256:,d} qubits for 1 run "
              f"(vs {v256:,d} for the plain monolithic QUBO).")
    out['single_run_collision'] = {m: 4 * V_min(m, 256)[0] for m in ('binary', 'wallace')}

    json.dump(out, open(os.path.join(_HERE, 'sampling.json'), 'w'), indent=1)
    print("\nwrote multirun/sampling.json")
