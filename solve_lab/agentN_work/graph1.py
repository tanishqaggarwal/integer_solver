"""Connected components of (a) eq-var bipartite graph, (b) eq-atom bipartite graph."""
import model, sys
from collections import Counter, defaultdict

d = model.get()
eq_terms = d['eq_terms']; atom_vars = d['atom_vars']; NA=len(atom_vars); NE=len(eq_terms)

class DSU:
    def __init__(s,n): s.p=list(range(n))
    def f(s,x):
        while s.p[x]!=x: s.p[x]=s.p[s.p[x]]; x=s.p[x]
        return x
    def u(s,a,b):
        a,b=s.f(a),s.f(b)
        if a!=b: s.p[a]=b

NV = 38748
# (a) variables joined when they co-occur in an equation
dsu = DSU(NV)
eqvars = []
for m,sq,tl in eq_terms:
    vs = set()
    for c,a in tl: vs.update(atom_vars[a])
    vs = sorted(vs); eqvars.append(vs)
    for x in vs[1:]: dsu.u(vs[0], x)
comp = defaultdict(list)
used = set()
for vs in eqvars: used.update(vs)
for x in used: comp[dsu.f(x)].append(x)
sizes = sorted((len(v) for v in comp.values()), reverse=True)
print('vars appearing in equations:', len(used), 'of', NV)
print('num components:', len(sizes), 'sizes top20:', sizes[:20])
print('size hist:', Counter(sizes).most_common(10))

# (b) atoms joined when they co-occur in an equation
dsu2 = DSU(NA)
for m,sq,tl in eq_terms:
    aa = [a for c,a in tl]
    for a in aa[1:]: dsu2.u(aa[0], a)
comp2 = defaultdict(list)
for a in range(NA): comp2[dsu2.f(a)].append(a)
s2 = sorted((len(v) for v in comp2.values()), reverse=True)
print('atom components:', len(s2), 'top20:', s2[:20])
print('atom size hist:', Counter(s2).most_common(10))

# equation sizes
print('eq varcount hist:', Counter(len(v) for v in eqvars).most_common(10))
print('max eq varcount:', max(len(v) for v in eqvars))
