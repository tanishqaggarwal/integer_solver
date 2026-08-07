#!/usr/bin/env python3
"""EXHAUSTIVE link (ii): for every one of the 256 conditional-pin booleans, prove (by mod-p affine
rigidity, not by numerical repair) that turning it on forces the selected wire pair to that boolean's
two pin constants mod p.  Checkpointed after every boolean and resumable."""
import sys, os, json, pickle, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cfg_rigid2 import build, p

pins = json.load(open(os.path.join(HERE, 'pins.json')))
sup = pickle.load(open(os.path.join(HERE, 'supp.pkl'), 'rb'))
A = [b for b in sup['7715'] if str(b) in pins]
B = [b for b in sup['34554'] if str(b) in pins]
OUT = os.path.join(HERE, 'sweep_ii.json')
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
COORD = {'A': (12186, 16742), 'B': (14853, 24908)}


def one(bj, tree):
    partner = 5090 if tree == 'A' else 22106
    R = build([bj, partner])
    rel = R['rel']; C = R['CONST']
    out = {'tree': tree, 'conflicts': R['conflicts']}
    vals = []
    for w in COORD[tree]:
        r, a, b = rel(w)
        vals.append(str((a * C[r] + b) % p) if r in C else None)
    out['forced'] = vals
    cs = {str(c % p) for _, c in pins[str(bj)]}
    out['pinconsts'] = sorted(cs)
    out['all_forced'] = None not in vals
    out['match'] = out['all_forced'] and set(vals) <= cs
    return out


todo = [(b, 'A') for b in A] + [(b, 'B') for b in B]
t0 = time.time()
for bj, tree in todo:
    if str(bj) in res:
        continue
    try:
        r = one(bj, tree)
    except Exception as ex:
        r = {'tree': tree, 'error': repr(ex)}
    res[str(bj)] = r
    json.dump(res, open(OUT, 'w'))
    n = len(res)
    good = sum(1 for x in res.values() if x.get('match'))
    conf = sum(1 for x in res.values() if x.get('conflicts'))
    print('%3d/%d bit %-6s tree %s match=%-5s conflicts=%s   running: match %d, conflict %d, t=%.0fs'
          % (n, len(todo), bj, tree, r.get('match'), r.get('conflicts'), good, conf, time.time() - t0),
          flush=True)
bad = [k for k, x in res.items() if not x.get('match')]
print('DONE  %d/%d booleans confirm link (ii);  DISAGREEING: %s' % (len(res) - len(bad), len(res), bad))
