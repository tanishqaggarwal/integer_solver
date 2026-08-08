#!/usr/bin/env python3
"""run_256.py -- sound presolve (propagation) on the real 256-bit modmul,
unconditional and conditioned on k pinned operand bits (the interval-split
setting).  Propagation is a necessary-consequence closure => sound at any size;
its soundness was verified exactly against the enumerated ceiling at p<=251."""
import os, sys, json, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common, persist
from collections import Counter

KW = dict(mult='schoolbook', leaf=32, red='naf', mode='wallace')


def build():
    t = time.time()
    m = common.build_mm(256, p=common.real_p256(), **KW)
    print(f"built 256-bit modmul: {m['Q'].n} vars in {time.time()-t:.1f}s")
    return m


def prop_count(Q, fixed):
    dom = persist.prop_fix(Q, fixed)
    return dom


def main():
    m = build()
    Q, A, B = m['Q'], m['A'], m['B']
    rng = random.Random(7)

    # unconditional
    dom0 = prop_count(Q, {})
    print(f"\nUNCONDITIONAL propagation: {len(dom0)} / {Q.n} fixed "
          f"({100*len(dom0)/Q.n:.3f}%)  kinds={dict(Counter(Q.kind[v] for v in dom0))}")

    # conditioned on k pinned low bits of operand A (random target), and also
    # the all-zero pattern (which cascades hardest), for a bracket.
    print("\nCONDITIONED propagation vs # pinned low bits of operand A:")
    print(f"  {'k':>4} {'zero-pin fixed':>15} {'rand-pin fixed':>15} {'per-bit(rand)':>14}")
    ks = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256]
    rows = []
    atgt = rng.getrandbits(256)
    for k in ks:
        pin_zero = {A.bits[t]: 0 for t in range(k)}
        pin_rand = {A.bits[t]: (atgt >> t) & 1 for t in range(k)}
        dz = prop_count(Q, pin_zero)
        dr = prop_count(Q, pin_rand)
        nz, nr = len(dz) - 0, len(dr)
        perbit = (nr - len(dom0)) / k if k else 0
        rows.append(dict(k=k, zero=nz, rand=nr))
        print(f"  {k:>4} {nz:>15,} {nr:>15,} {perbit:>14.2f}")

    # both operands fully pinned => everything must be determined (sanity)
    allA = {A.bits[t]: (atgt >> t) & 1 for t in range(256)}
    btgt = rng.getrandbits(256)
    allAB = dict(allA); allAB.update({B.bits[t]: (btgt >> t) & 1 for t in range(256)})
    dAB = prop_count(Q, allAB)
    print(f"\n  A and B BOTH fully pinned: {len(dAB):,} / {Q.n:,} fixed "
          f"({100*len(dAB)/Q.n:.1f}%)")

    # the 'product known' factoring setting for contrast: pin C, operands free
    C = m['C']
    ctgt = rng.getrandbits(256)
    pinC = {C.bits[t]: (ctgt >> t) & 1 for t in range(256)}
    dC = prop_count(Q, pinC)
    print(f"  C (product) fully pinned, A,B free: {len(dC):,} / {Q.n:,} fixed "
          f"({100*len(dC)/Q.n:.3f}%)")

    with open(os.path.join(os.path.dirname(__file__), 'persist_256.json'), 'w') as f:
        json.dump(dict(n=Q.n, uncond=len(dom0), rows=rows,
                       both_pinned=len(dAB), c_pinned=len(dC)), f, indent=2)


if __name__ == '__main__':
    sys.setrecursionlimit(1000000)
    main()
