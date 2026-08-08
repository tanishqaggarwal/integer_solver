#!/usr/bin/env python3
"""rank.py -- assemble measured marginal window costs into whole-instance totals.

Same accounting as the existing report.py:  total = ceil(256 / w) * (cost of one
marginal window).  The marginal window already excludes the incoming x1,y1
registers, so the windows tile exactly; the leading accumulator and the two final
target congruences add ~1.5k qubits, under 0.05% of any total here, and are
omitted on both sides so the numbers line up with report.py's 9,061,804.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = 9061804          # report.py's headline: binary, w=9, random 256-bit p

LABEL = {
    'base': 'A0  baseline (weighted-sum table, schoolbook, chunk 16)',
    'base_c': 'A1  baseline + chunk retune                (same arithmetic)',
    'mux': 'B   + one-hot MUX table look-up',
    'muxtree': 'C   + AND-tree one-hot (no cardinality penalty)',
    'muxtreekara': 'D   + Karatsuba x4',
    'full': 'E   + signed digits',
    'fullc64': 'F   + chunk 64',
    'unsc64': 'F-  E/F without signed digits',
    'semaev': 'X   x-only Semaev S_3 chain (per S_3 step)',
    'semaevc64': 'X   x-only Semaev S_3 chain, chunk 64',
    'basekara': 'D-  baseline table + Karatsuba x4 only',
    'best': 'G   + Toom-3 x2 + pseudo-Mersenne reduction (REAL p)',
    'bestc128': 'G   + Toom-3 x2 + pseudo-Mersenne, chunk 128 (REAL p)',
    'realbase': 'A0r baseline arithmetic at the REAL p',
    'realmux': 'B r  MUX+tree at the REAL p',
}


def load():
    W = json.load(open(os.path.join(HERE, 'win256.json')))
    fam = {}
    for k, v in W.items():
        if '_w' not in k:
            continue
        f, w = k.rsplit('_w', 1)
        w = int(w.split('_')[0])
        suffix = k.rsplit('_w', 1)[1]
        if '_' in suffix:                       # e.g. full_w12_c32
            f = f + '_' + suffix.split('_', 1)[1]
        fam.setdefault(f, []).append((w, v))
    return fam


def main():
    fam = load()
    rows = []
    print("=" * 104)
    print("MARGINAL COST OF ONE COMB WINDOW AT s = 256, AND THE WHOLE-INSTANCE TOTAL")
    print("=" * 104)
    print(f"{'family':<26}{'w':>4}{'windows':>9}{'qubits/window':>15}"
          f"{'TOTAL qubits':>15}{'TOTAL couplers':>17}{'|J|':>7}{'vs base':>9}")
    for f in sorted(fam):
        best = None
        for w, v in sorted(fam[f]):
            M = math.ceil(256 / w)
            tv, tc = M * v['vars'], M * v['couplers']
            print(f"{f:<26}{w:>4}{M:>9}{v['vars']:>15,d}{tv:>15,d}{tc:>17,d}"
                  f"{'2^%d' % v['dynamic_range_bits']:>7}{BASE / tv:>8.2f}x")
            if best is None or tv < best[1]:
                best = (w, tv, tc, v['dynamic_range_bits'], v['vars'])
        rows.append((f, best))
        print()
    print("=" * 104)
    print("BEST w PER FAMILY, RANKED")
    print("=" * 104)
    print(f"{'family':<50}{'w':>4}{'TOTAL qubits':>15}{'couplers':>16}"
          f"{'|J|':>7}{'vs 9.06e6':>11}")
    for f, b in sorted(rows, key=lambda r: r[1][1]):
        print(f"{LABEL.get(f, f):<50}{b[0]:>4}{b[1]:>15,d}{b[2]:>16,d}"
              f"{'2^%d' % b[3]:>7}{BASE / b[1]:>10.2f}x")
    print()
    best = min(b[1] for f, b in rows if not f.startswith('semaev'))
    for n, N in (("D-Wave Advantage  (5,760 qubits)", 5760),
                 ("D-Wave Advantage2 (4,400 qubits)", 4400)):
        print(f"{n}: best exact encoding is {best / N:,.0f}x too big")


if __name__ == '__main__':
    main()
