#!/usr/bin/env python3
"""compose.py -- window-width co-optimization of the FULL 256-bit comb ladder.

Composes the measured per-window cost (real p, mm.py live) into the whole ladder:

    ladder(w) = ceil(256/w) * window(w)

where window(w) is the lab's own marginal comb window (squeeze/window.py): two
one-hot table look-ups, four linear words, THREE modular multiplications
(lam*d general, lam*lam SQUARING, lam*(x1-x3) general) and two d!=0 gadgets,
with the carried-in x1,y1 (2*s bits) subtracted so the units compose.

Physical qubits from the whole-ladder max clique K:
    physical = logical * max(1, K/12)   on Pegasus  (Advantage, P16)
             = logical * max(1, K/16)   on Zephyr   (Advantage2, Z6)
A clique K_c is embedded with chains of length ~c/12 (Pegasus) / ~c/16 (Zephyr);
for K<=12/16 the chains are length 1, so physical == logical (no overhead).
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.normpath(os.path.join(HERE, '..', '..', 'squeeze'))
sys.path.insert(0, SQ)
from window import window                                        # noqa: E402

HW_PEG, HW_ZEP = 5760, 4400


def phys(logical, K, div):
    return logical * max(1.0, K / div)


def main():
    cache = os.path.join(HERE, 'windows_full.json')
    res = json.load(open(cache)) if os.path.exists(cache) else {}
    rows = []
    for w in range(4, 13):
        key = f"w{w}"
        if key not in res:
            res[key] = window(w, mode='wallace', mult='karatsuba', leaf=24,
                              red='naf', naf_consts=True)
            json.dump(res, open(cache, 'w'), indent=1, sort_keys=True)
        r = res[key]
        M = math.ceil(256 / w)
        log = M * r['vars']
        coup = M * r['couplers']
        rows.append(dict(w=w, M=M, win=r['vars'], K=r['clique'],
                         jbits=r['jbits'], win_coup=r['couplers'],
                         logical=log, couplers=coup))
    print("PER-WINDOW cost (real p, s=256, karatsuba(24)/naf/wallace)")
    print(f"  {'w':>3} {'D=2^w':>7} {'window vars':>12} {'K':>3} {'|J|':>5} "
          f"{'window couplers':>16}")
    for r in rows:
        print(f"  {r['w']:3d} {1<<r['w']:7d} {r['win']:12,d} {r['K']:3d} "
              f"2^{r['jbits']:<3d} {r['win_coup']:16,d}")

    print("\nFULL 256-bit LADDER  = ceil(256/w) * window(w)")
    hdr = (f"  {'w':>3} {'windows':>7} {'LOGICAL':>14} {'K':>3} {'|J|':>5} "
           f"{'couplers':>15} {'phys(Peg)':>14} {'phys(Zep)':>14}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    best = None
    for r in rows:
        pp = phys(r['logical'], r['K'], 12)
        pz = phys(r['logical'], r['K'], 16)
        print(f"  {r['w']:3d} {r['M']:7d} {r['logical']:14,d} {r['K']:3d} "
              f"2^{r['jbits']:<3d} {r['couplers']:15,d} {pp:14,.0f} {pz:14,.0f}")
        if best is None or r['logical'] < best['logical']:
            best = dict(r, phys_peg=pp, phys_zep=pz)
    print("  " + "-" * (len(hdr) - 2))
    print(f"\nSMALLEST FULL LADDER: w={best['w']} ({best['M']} windows)")
    print(f"  logical      = {best['logical']:,d}")
    print(f"  max clique K = {best['K']}")
    print(f"  |J|          = 2^{best['jbits']}")
    print(f"  couplers     = {best['couplers']:,d}")
    print(f"  physical Pegasus (Advantage,  5760 q) = {best['phys_peg']:,.0f}"
          f"  = {best['phys_peg']/HW_PEG:,.0f}x chip")
    print(f"  physical Zephyr  (Advantage2, 4400 q) = {best['phys_zep']:,.0f}"
          f"  = {best['phys_zep']/HW_ZEP:,.0f}x chip")
    json.dump({'rows': rows, 'best': best},
              open(os.path.join(HERE, 'ladder.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
