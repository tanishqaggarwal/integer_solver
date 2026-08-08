#!/usr/bin/env python3
"""crossover.py -- where recursion starts to pay, measured not estimated.

For each width s, build ONE modular multiplication mod a pseudo-Mersenne prime
of that width, schoolbook vs Karatsuba vs Toom-3, and print the measured
(logical, clique, |J|) triple.  The crossover is the smallest s at which the
recursive form is cheaper than schoolbook.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measure import modmul                                       # noqa: E402

# pseudo-Mersenne primes 2^k - c of each width
PM = {8: 2**8 - 5, 12: 2**12 - 3, 16: 2**16 - 15, 24: 2**24 - 3,
      32: 2**32 - 5, 48: 2**48 - 59, 64: 2**64 - 59, 96: 2**96 - 17,
      128: 2**128 - 159, 192: 2**192 - 237, 256: 2**256 - 2**32 - 977}


def best(mult, s, p, mode, leaves):
    out = None
    for lf in leaves:
        if lf >= s:
            continue
        r = modmul(mult=mult, leaf=lf, red='naf', mode=mode, s=s, p=p)
        r['leaf'] = lf
        if out is None or r['vars'] < out['vars']:
            out = r
    return out


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'wallace'
    leaves = (4, 8, 12, 16, 24, 32, 64)
    print(f"one modular multiplication, mode={mode}, reduction=naf")
    print(f"{'s':>5} {'schoolbook':>12} {'karatsuba':>12} {'leaf':>5} {'ratio':>7} "
          f"{'toom3':>12} {'leaf':>5} {'ratio':>7}")
    for s in sorted(PM):
        p = PM[s]
        sb = modmul(mult='schoolbook', red='naf', mode=mode, s=s, p=p)
        ka = best('karatsuba', s, p, mode, leaves)
        tm = best('toom3', s, p, mode, leaves)
        row = f"{s:5d} {sb['vars']:12,d}"
        for r in (ka, tm):
            if r is None:
                row += f" {'-':>12} {'-':>5} {'-':>7}"
            else:
                row += f" {r['vars']:12,d} {r['leaf']:5d} {r['vars']/sb['vars']:7.3f}"
        print(row)
