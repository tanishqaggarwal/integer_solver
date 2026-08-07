#!/usr/bin/env python3
"""|S| <= 4 probe for all 12 orbit targets (exact: truncation gives false positives only)."""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'agentX_work', 'pylib'))
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'ydata.json')))
L = [(int(x), int(y)) for x, y in d['ladder']]
orb = json.load(open(os.path.join(HERE, 'yorbit.json')))
k = np.memmap(os.path.join(HERE, '..', 'agentX_work', 'tbl4s.bin'), dtype=np.uint64, mode='r')
def has(q):
    q = np.uint64(q); i = int(np.searchsorted(k, q))
    return i < len(k) and int(k[i]) == int(q)
print('%-10s  %-8s  %-8s' % ('target', '|S|=1?', '|S|2..4?'))
for nm, P in sorted(orb.items()):
    Px = int(P[0]); Py = int(P[1])
    one = [i for i in range(256) if L[i] == (Px, Py)]
    print('%-10s  %-8s  %-8s' % (nm, one if one else 'no', 'HIT' if has(Px & 0xFFFFFFFFFFFFFFFF) else 'no'))
