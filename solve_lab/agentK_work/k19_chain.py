#!/usr/bin/env python3
"""K19: (a) order the 256 leaves into the doubling chain, (b) split them by which root
half their wire feeds."""
import sys, os, json, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from cascadep import CascadeP, NV, P
from circ2 import vars_of

D = FD.points()
L = D['leaves']
pts = [(int(l['X']), int(l['Y'])) for l in L]
sels = [l['sel'] for l in L]
S = {p: i for i, p in enumerate(pts)}
dbl = {}
for i, p in enumerate(pts):
    q = FD.add(p, p)
    if q in S: dbl[i] = S[q]
isdouble = set(dbl.values())
starts = [i for i in range(256) if i not in isdouble]
print('chain starts (not the double of any leaf):', starts)
ends = [i for i in range(256) if i not in dbl]
print('chain ends (double not in the set):', ends)
assert len(starts) == 1 and len(ends) == 1
chain = [starts[0]]
while chain[-1] in dbl: chain.append(dbl[chain[-1]])
print('chain length:', len(chain))
exp = {chain[k]: k for k in range(len(chain))}   # leaf index -> exponent
assert len(exp) == 256

# (b) which root half does each leaf feed?
C = CascadeP()
defnode = {}
for a in C.E.order:
    c = C.E.cls[a]; defnode[c[1]] = c[2]
uses = collections.defaultdict(list)
for w, n in defnode.items():
    for z in vars_of(n):
        if z != w: uses[z].append(w)

def reaches(start, goals, cap=400000):
    seen = set([start]); st = [start]; hit = set()
    while st:
        u = st.pop()
        if u in goals: hit.add(u); continue
        for w in uses.get(u, ()):
            if w not in seen:
                seen.add(w); st.append(w)
        if len(seen) > cap: break
    return hit

GOALS = {23927: 'A', 19083: 'A', 1308: 'B', 17601: 'B'}
side = {}
for i, l in enumerate(L):
    h = reaches(l['wx'], set(GOALS))
    tags = set(GOALS[g] for g in h)
    side[i] = ''.join(sorted(tags))
cnt = collections.Counter(side.values())
print('side tags:', cnt)
A = [i for i in range(256) if side[i] == 'A']
B = [i for i in range(256) if side[i] == 'B']
print('A-side leaves', len(A), 'B-side leaves', len(B))
print('A exponents (sorted):', sorted(exp[i] for i in A)[:20], '...')
print('B exponents (sorted):', sorted(exp[i] for i in B))
json.dump({'exp': {str(i): exp[i] for i in range(256)},
           'sel': {str(i): sels[i] for i in range(256)},
           'side': {str(i): side[i] for i in range(256)},
           'A_exp': sorted(exp[i] for i in A), 'B_exp': sorted(exp[i] for i in B)},
          open(K + '/chain.json', 'w'))
