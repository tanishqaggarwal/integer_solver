#!/usr/bin/env python3
"""Reduce the atom system mod p = 2^256-2^32-977 and propagate from the pins.

Over F_p every gate is invertible, so propagation is far stronger than over Z.
"""
import pickle, os, collections, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
P = 2**256 - 2**32 - 977
NV = 38748


def load():
    D = pickle.load(open(os.path.join(HERE, 'atoms.pkl'), 'rb'))
    polys = pickle.load(open(os.path.join(HERE, 'polys.pkl'), 'rb'))
    mp = []
    for q in polys:
        r = {}
        for m, c in q.items():
            c %= P
            if c:
                r[m] = c
        mp.append(r)
    return D, polys, mp


if __name__ == '__main__':
    D, polys, mp = load()
    tri = sum(1 for r in mp if not r)
    print("atoms that vanish identically mod p:", tri)
    import re
    sh = collections.Counter()
    for i, r in enumerate(mp):
        vs = set()
        for m in r:
            vs |= set(m)
        deg = max((len(m) for m in r), default=0)
        sh[(deg, len(r), len(vs))] += 1
    print("mod-p (deg,nterms,nvars) histogram:")
    for k, v in sorted(sh.items()):
        print("   ", k, v)
