#!/usr/bin/env python3
"""Scan the 256 conditional-pin selectors: force each ON (and OFF) and measure
the mod-p state: conflicts, A = X35389, B = X6671, coordinate values."""
import pickle, os, sys, re, time, json, collections
from boolscore import Fast
from fp import P
HERE = os.path.dirname(os.path.abspath(__file__))

F = Fast(); M = F.M
sel = collections.defaultdict(list)
for i, s in enumerate(M.src):
    m = re.match(r'^X(\d+) \* \(X(\d+) - (\d+)\)', s)
    if m and int(m.group(3)) > 2**200:
        sel[int(m.group(1))].append((int(m.group(2)), int(m.group(3))))
selectors = sorted(sel)
print("selectors:", len(selectors), flush=True)
COORD = {'x1': 12186, 'y1': 16742, 'x2': 14853, 'y2': 24908,
         'x3': 22162, 'y3': 30213, 'A': 35389, 'B': 6671}
lo = int(sys.argv[1]) if len(sys.argv) > 1 else 0
hi = int(sys.argv[2]) if len(sys.argv) > 2 else len(selectors)
mode = sys.argv[3] if len(sys.argv) > 3 else 'on'
out = {}
t0 = time.time()
path = os.path.join(HERE, f'sel_{mode}_{lo}_{hi}.json')
for k, b in enumerate(selectors[lo:hi]):
    tgtval = 1 if mode == 'on' else 0

    def pol(v, roots, tgt=b, tv=tgtval):
        if v == tgt:
            return tv if tv in roots else roots[0]
        x = F.witp[v]
        return x if x in roots else roots[0]

    val, conf, dd = F.run(pol)
    rec = {'nconf': len(conf), 'conf': conf[:16]}
    for nm, v in COORD.items():
        rec[nm] = val[v]
    out[b] = rec
    if k % 10 == 0:
        json.dump({str(a): c for a, c in out.items()}, open(path, 'w'))
        print(f"{lo+k}/{hi} sel X{b} nconf={len(conf)} A={val[35389]} t={time.time()-t0:.0f}s",
              flush=True)
json.dump({str(a): c for a, c in out.items()}, open(path, 'w'))
print("done", time.time() - t0)
