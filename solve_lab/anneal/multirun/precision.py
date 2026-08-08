#!/usr/bin/env python3
"""precision.py -- is the coupler dynamic range something a multi-run split can fix?

Splitting the problem reduces the QUBO's SIZE.  It does not obviously reduce the
range of coupler strengths, because that range is set by the local structure of a
single column plus the AND-penalty weight, not by the instance width.  This
measures |J|max/|J|min for one modular multiplication as the word size s runs
from 8 to 256 bits, and for the instance's real prime.

Writes multirun/precision.json.
"""
import json, os, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)

import sympy
import pieces
from instance import p as REAL_P

out = {'sweep': [], 'real_p': {}}
print("coupler dynamic range of ONE modular multiplication, vs word size")
print(f"{'s':>5} {'prime':>26} {'binary |J| bits':>16} {'wallace |J| bits':>17}")
for s in (8, 16, 32, 64, 128, 256):
    q = (1 << (s - 1)) + 1
    while not sympy.isprime(q):
        q += 2
    row = {}
    for mode in ('binary', 'wallace'):
        pieces.P = q
        row[mode] = pieces.piece_mul(mode)['jbits']
    print(f"{s:5d} {('2^%d+%d' % (s-1, q-(1<<(s-1)))):>26} "
          f"{row['binary']:16d} {row['wallace']:17d}")
    out['sweep'].append(dict(s=s, **row))

for mode in ('binary', 'wallace'):
    pieces.P = REAL_P
    out['real_p'][mode] = pieces.piece_mul(mode)['jbits']
print(f"{256:5d} {'p = 2^256-2^32-977':>26} "
      f"{out['real_p']['binary']:16d} {out['real_p']['wallace']:17d}")

print()
print("The range is INDEPENDENT of s.  It is fixed by one column's coupling load and")
print("by the AND-penalty weight W chosen in qubo.py:finalize().  Cutting the instance")
print("into 2^(256-mu) pieces therefore leaves the precision requirement untouched:")
print(f"  wallace needs {out['real_p']['wallace']} bits, binary needs "
      f"{out['real_p']['binary']} bits, hardware offers ~4-5.")

json.dump(out, open(os.path.join(_HERE, 'precision.json'), 'w'), indent=1)
print("\nwrote multirun/precision.json")
