#!/usr/bin/env python3
"""mincost_modmul.py -- push the ~120k-qubit 256-bit modular multiply smaller.
Baseline (squeeze): karatsuba leaf-24/naf/wallace = 99,298 (clique 5, |J| 2^6).
Carries are 2/3 of the cost -> the summation network is the lever.
Plan sweep: toom3(128)>kara(24) = 96,809."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'squeeze'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from measure import modmul
PLANS = [('kara-24 (squeeze best)', dict(mult='karatsuba', leaf=24)),
         ('toom3(128)>kara(24)', dict(mult=[('toom3',128),('karatsuba',24)], leaf=24)),
         ('toom3(128)>kara(16)', dict(mult=[('toom3',128),('karatsuba',16)], leaf=16))]
if __name__ == '__main__':
    print(f"{'plan':26} {'vars':>7} {'clique':>6} {'|J|':>5} {'AND':>7} {'carry':>7}")
    for lab, kw in PLANS:
        st = modmul(red='naf', mode='wallace', **kw)
        print(f"{lab:26} {st['vars']:7d} {st['clique']:6d} 2^{st['jbits']:<3d} {st['ands']:7d} {st['carries']:7d}")
