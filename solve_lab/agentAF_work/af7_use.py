#!/usr/bin/env python3
"""agent AF, step 7: where the P*u slack wires are USED -> the actual lift conditions."""
import sys, os, pickle
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import (atoms, defs, defc, val, lift, Pval, find, pp, expand, varsof, shape_of)

real = [t for t in lift if t[4] == Pval]
Rset = {}
for (aid, Rw, mv, uv, M) in real:
    Rset.setdefault(find(Rw), []).append((aid, find(uv)))
print('distinct slack wires R (classes): %d  from %d atoms' % (len(Rset), len(real)))
print('R wires used by >1 lift atom:', sum(1 for k, v in Rset.items() if len(v) > 1))

# index: variable class -> atoms mentioning it
v2a = defaultdict(list)
for aid, a in enumerate(atoms):
    for v in varsof(a, set()):
        v2a[find(v)].append(aid)

pat = Counter()
uses = {}
for r in Rset:
    lst = [aid for aid in v2a[r] if aid not in {x[0] for x in Rset[r]}]
    uses[r] = lst
    pat[len(lst)] += 1
print('non-defining uses per slack wire:', dict(sorted(pat.items())))

sh = Counter()
for r, lst in uses.items():
    for aid in lst:
        sh[shape_of(atoms[aid])] += 1
print('\nshapes of the USE atoms:')
for s, k in sh.most_common(25):
    print('  %6d  %s' % (k, s))

print('\nexamples:')
n = 0
for r, lst in uses.items():
    if n >= 8:
        break
    for aid in lst:
        print('  R=x%d  use atom %d : %s' % (r, aid, pp(atoms[aid])[:220]))
    n += 1
