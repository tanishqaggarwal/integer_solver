#!/usr/bin/env python3
"""agent AF, step 16: incidence of the 766 off-pin free wires, and the congruence residual chain."""
import sys, os, pickle
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
conds = C['conds']; info = M['info']

v2a = defaultdict(list)
for aid, a in enumerate(atoms):
    for v in varsof(a, set()):
        v2a[find(v)].append(aid)

opw = [(i, conds[i][1], info[i]['other'], info[i]['gate'])
       for i in range(len(conds)) if info[i]['cls'] == 'offpin']
print('off-pin conditions: %d ; distinct free wires: %d' % (len(opw), len(set(t[2] for t in opw))))
inc = Counter(len(v2a[w]) for (_, _, w, _) in opw)
print('#atoms mentioning an off-pin wire:', dict(sorted(inc.items())))

sh = Counter()
for (_, _, w, g) in opw:
    for aid in v2a[w]:
        sh[shape_of(atoms[aid])] += 1
print('\nshapes of atoms touching off-pin wires:')
for s, k in sh.most_common(14):
    print('   %6d  %s' % (k, s))

# For each off-pin wire, are ALL its other occurrences multiplied by the SAME gate L?
def raw(n):
    if n[0] == 'c':
        return str(n[1]) if abs(n[1]) < 10**14 else ('P' if n[1] == Pval else 'BIG')
    if n[0] == 'v':
        return 'x%d' % find(n[1])
    if n[0] == 'neg':
        return '-' + raw(n[1])
    return '(%s %s %s)' % (raw(n[1]), n[0], raw(n[2]))
print('\nall atoms touching one off-pin wire (3 examples):')
for (i, c, w, g) in opw[:3]:
    print('  wire x%d   (block gate x%d, c=%d)' % (w, g, c))
    for aid in v2a[w]:
        print('      ', raw(atoms[aid])[:200])
