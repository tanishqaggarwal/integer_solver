#!/usr/bin/env python3
"""atom.py -- freshly measure the best modmul ATOM at s=256, importing mm LIVE.

Does NOT read squeeze/results.json (which is cached against an older mm.py).
Sweeps the Karatsuba leaf so we report the CURRENT best general modmul and the
current best squaring, in the hardware currency (logical, K, |J|, couplers).
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.normpath(os.path.join(HERE, '..', '..', 'squeeze'))
sys.path.insert(0, SQ)
from mmqb import MMQB                                             # noqa: E402
from mm import build_modmul                                      # noqa: E402

P = 2 ** 256 - 2 ** 32 - 977
S = 256


def atom(square=False, mult='karatsuba', leaf=24, red='naf', mode='wallace'):
    Q = MMQB(chunk=16, mode=mode)
    A = Q.mkword('A', S, lambda wv: wv['_a'])
    B = A if square else Q.mkword('B', S, lambda wv: wv['_b'])
    C = Q.mkword('C', S, lambda wv: wv['_c'])
    base = Q.n
    t0 = time.time()
    build_modmul(Q, P, A, B, C, mult=mult, leaf=leaf, red=red)
    Q.finalize()
    st = Q.stats()
    return dict(vars=Q.n - base, clique=st['max_clique'],
                jbits=st['dynamic_range_bits'], couplers=st['couplers'],
                ands=st['and_vars'], secs=round(time.time() - t0, 1))


if __name__ == '__main__':
    print("FRESH ATOM MEASUREMENT (s=256, real p=2^256-2^32-977, mm.py imported live)")
    print(f"{'atom':>40} {'logical':>10} {'K':>4} {'|J|':>6} {'couplers':>12} "
          f"{'ANDs':>8} {'s':>6}")
    print("-" * 92)
    best = {}
    for sq in (False, True):
        label = 'SQUARING' if sq else 'GENERAL '
        for leaf in (16, 18, 20, 24, 28, 32):
            r = atom(square=sq, leaf=leaf)
            key = f"{label} karatsuba(leaf={leaf})/naf/wallace"
            print(f"{key:>40} {r['vars']:10,d} {r['clique']:4d} 2^{r['jbits']:<3d} "
                  f"{r['couplers']:12,d} {r['ands']:8,d} {r['secs']:6}")
            if label not in best or r['vars'] < best[label][1]['vars']:
                best[label] = (leaf, r)
    print("-" * 92)
    import json
    out = {}
    for label, (leaf, r) in best.items():
        print(f"BEST {label}: leaf={leaf}  {r['vars']:,d} logical, K={r['clique']}, "
              f"|J|=2^{r['jbits']}, {r['couplers']:,d} couplers")
        out[label.strip()] = dict(leaf=leaf, **r)
    json.dump(out, open(os.path.join(HERE, 'atom.json'), 'w'), indent=1)
