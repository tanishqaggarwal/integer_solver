#!/usr/bin/env python3
"""Measure the effect of flipping each boolean decision on the mod-p state."""
import pickle, os, sys, time, json
from boolscore import Fast
from fp import P
HERE = os.path.dirname(os.path.abspath(__file__))

F = Fast()
base = pickle.load(open(os.path.join(HERE, 'bool_wit.pkl'), 'rb'))
decs = base['dec']
lo = int(sys.argv[1]); hi = min(int(sys.argv[2]), len(decs))
out = {}
t0 = time.time()
outpath = os.path.join(HERE, f'flip_{lo}_{hi}.json')
for k, u in enumerate(decs[lo:hi]):
    w = F.witp[u]
    if w not in (0, 1):
        continue
    other = 1 - w

    def pol(v, roots, tgt=u, o=other):
        if v == tgt:
            return o if o in roots else roots[0]
        x = F.witp[v]
        return x if x in roots else roots[0]

    val, conf, dd = F.run(pol)
    out[u] = {'nconf': len(conf), 'conf': conf[:16],
              'A': val[35389], 'B': val[6671],
              'known': sum(1 for x in val if x is not None)}
    if k % 20 == 0:
        json.dump({str(a): b for a, b in out.items()}, open(outpath, 'w'))
        print(f"{lo+k}/{hi} X{u} nconf={len(conf)} t={time.time()-t0:.0f}s", flush=True)
json.dump({str(a): b for a, b in out.items()}, open(outpath, 'w'))
print("done", lo, hi, time.time() - t0)
