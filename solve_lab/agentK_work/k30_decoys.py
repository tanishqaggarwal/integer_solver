#!/usr/bin/env python3
"""K30: are the 900 non-leaf free booleans real knobs or dead ends?

My §1 says the only boolean inputs are the 256 leaf selectors.  The fixed classifier
(k25_class.py) finds 1156 free booleans, i.e. 900 besides the leaf selectors.  If any of
those 900 can influence a leaf wire, a slot wire or the root, then §1's knob set is wrong
and the search space is bigger than I claimed.  This decides it by forward cone."""
import sys, os, json, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P
from circ2 import vars_of

C = CascadeP()
vc = json.load(open(K + '/varclass2.json'))
other = vc['otherbools']
print('non-leaf free booleans:', len(other))

D = json.load(open(K + '/points.json'))
LEAFW = set()
for l in D['leaves']: LEAFW.add(l['wx']); LEAFW.add(l['wy'])
ROOTV = {12186, 16742, 14853, 24908, 22162, 30213, 24468, 18956, 13682, 37892, 15298}
TARGETS = LEAFW | ROOTV

# undirected atom graph: two variables are adjacent if they share an atom.  This is the
# MOST generous notion of influence - if a knob cannot reach a target even here, it is dead.
adj = collections.defaultdict(set)
for i, vs in enumerate(C.avars):
    if len(vs) > 12: continue           # skip the few very wide atoms
    for a in vs:
        for b in vs:
            if a != b: adj[a].add(b)

reach_sizes = []
live = []
for u in other:
    seen = {u}; st = [u]; hit = False
    while st:
        x = st.pop()
        if x in TARGETS: hit = True; break
        for y in adj.get(x, ()):
            if y not in seen:
                seen.add(y); st.append(y)
        if len(seen) > 60: break         # anything this big is not a dead end
    reach_sizes.append(len(seen))
    if hit or len(seen) > 60: live.append((u, len(seen), hit))

h = collections.Counter(reach_sizes)
print('cone-size histogram (undirected, capped at 60):', sorted(h.items())[:15])
print('knobs whose cone reaches a leaf wire / root var, or exceeds 60 vars:', len(live))
for u, s, hit in live[:15]:
    print('   x%d cone=%d reached_target=%s  atoms:' % (u, s, hit))
    for i in C.var2atoms[u][:4]: print('       ', C.names[i][:95])
json.dump({'live': [[u, s, hit] for u, s, hit in live]}, open(K + '/decoys.json', 'w'))
