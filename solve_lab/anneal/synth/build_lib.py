#!/usr/bin/env python3
"""build_lib.py -- generate a library of synthetic planted-key instances, cache to JSON."""
import sys, json, time, os
sys.path.insert(0, '.')
from synth.gen import make, Curve, Instance

LIB = 'synth/lib.json'

def load():
    if not os.path.exists(LIB): return {}
    return json.load(open(LIB))

def save(d): json.dump(d, open(LIB, 'w'), indent=1)

if __name__ == '__main__':
    sizes = [int(x) for x in sys.argv[1:]] or [8, 12, 16, 20, 24, 28, 32]
    d = load()
    for b in sizes:
        if str(b) in d: continue
        t = time.time(); inst = make(b, seed=3)
        d[str(b)] = dict(p=inst.curve.p, B=inst.curve.B, Gx=inst.G[0], Gy=inst.G[1],
                         n=inst.n, k=inst.k, Tx=inst.T[0], Ty=inst.T[1], bits=b)
        save(d)
        print(f"cached {b}-bit ({time.time()-t:.1f}s): k={inst.k}", flush=True)
    print("done")

def get(b):
    r = load()[str(b)]
    c = Curve(r['p'], r['B'])
    return Instance(c, (r['Gx'], r['Gy']), r['n'], r['k'], (r['Tx'], r['Ty']), r['bits'])
