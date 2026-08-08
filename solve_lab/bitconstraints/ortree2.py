#!/usr/bin/env python3
"""Pass 13: reconstruct the OR-tree over the 256 selector bits.

The OR gate is not a single atom.  It is a THREE-atom gadget:
    sum   :  s - u - v        ->  s = u + v
    prod  :  u*v - p          ->  p = u*v
    diff  :  s - p - o        ->  o = s - p  =  u + v - u*v  =  u OR v
(plus copy atoms  o - o'  that rename the wire).
We detect all such gadgets, build the OR DAG upward from the 256 selectors,
find its root, and check whether the root is pinned to 1 -- that would be the
cardinality constraint  OR(b_0..b_255) = 1, i.e. "not all selectors zero".
"""
import os, json
from collections import defaultdict, Counter
from scan import load, support, degree, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain, BITS = load()
atoms = D['atoms']
atom_eqs = defaultdict(list)
for i, terms in enumerate(D['eq_terms']):
    for c, aid in terms:
        atom_eqs[aid].append(i)

# ---- primitive gate extraction ----
lin3 = []          # (atom, w, u, v)  meaning  w = u + v
prod = defaultdict(list)   # frozenset{u,v} -> [(atom, p), ...]
copy_edges = defaultdict(list)
pin1 = {}
for i, a in enumerate(atoms):
    d = dict(a); s = support(a); dg = degree(a)
    if len(a) == 2 and len(s) == 1 and dg == 1 and d.get(()) == 1 and d.get((next(iter(s)),)) == -1:
        pin1.setdefault(next(iter(s)), i)
    elif len(a) == 2 and len(s) == 2 and dg == 1 and set(d.values()) == {1, -1}:
        u = [m[0] for m in d if d[m] == 1][0]
        v = [m[0] for m in d if d[m] == -1][0]
        copy_edges[u].append((v, i)); copy_edges[v].append((u, i))
    elif len(a) == 3 and len(s) == 3 and dg == 1 and all(abs(c) == 1 for c in d.values()):
        pos = [m[0] for m in d if d[m] == 1]
        neg = [m[0] for m in d if d[m] == -1]
        if len(pos) == 2:            #  pos0 + pos1 - neg0 = 0  =>  neg0 = pos0+pos1
            lin3.append((i, neg[0], pos[0], pos[1]))
        else:                        #  pos0 - neg0 - neg1 = 0  =>  pos0 = neg0+neg1
            lin3.append((i, pos[0], neg[0], neg[1]))
    elif len(a) == 2 and len(s) == 3 and dg == 2:
        q = [m for m in d if len(m) == 2]; l = [m for m in d if len(m) == 1]
        if len(q) == 1 and len(l) == 1 and d[q[0]] == -d[l[0]] and q[0][0] != q[0][1]:
            prod[frozenset(q[0])].append((i, l[0][0]))

print(f"sum atoms  w = u+v : {len(lin3)}")
print(f"prod atoms p = u*v : {len(prod)}")
print(f"copy atoms         : {sum(len(v) for v in copy_edges.values())//2}")
print(f"unit pins  x = 1   : {len(pin1)}")

sum_by_pair = defaultdict(list)      # {u,v} -> [(atom, w)]
sum_by_operands = defaultdict(list)  # ordered lookup: (w) -> list of (u,v,atom)
for i, w, u, v in lin3:
    sum_by_pair[frozenset((u, v))].append((i, w))
    sum_by_operands[w].append((u, v, i))

# ---- OR gadgets ----
ors = []                              # (o, u, v, [atom ids])
for pair, plist in prod.items():
    if len(pair) != 2:
        continue
    for (pa, p) in plist:
        for (sa, s) in sum_by_pair.get(pair, []):
            # look for  o = s - p, i.e. a sum atom  s = o + p
            for (x, y, da) in sum_by_operands.get(s, []):
                if x == p and y != p:
                    ors.append((y, tuple(pair), [sa, pa, da]))
                elif y == p and x != p:
                    ors.append((x, tuple(pair), [sa, pa, da]))
print(f"\nOR gadgets found (o = u OR v): {len(ors)}")

or_children = {}
for o, (u, v), prov in ors:
    or_children.setdefault(o, (u, v, prov))

# ---- alias classes through copy atoms ----
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb
for u, lst in copy_edges.items():
    for v, i in lst:
        union(u, v)

# ---- climb the OR DAG from the 256 selectors ----
# leaf set of each OR wire
leafset = {}
for b in chain:
    leafset[find(b)] = frozenset([b])
# constant-zero wires act as neutral OR leaves
zero_atoms = {next(iter(support(a))) for a in atoms
              if len(a) == 1 and len(support(a)) == 1 and degree(a) == 1}
for z in zero_atoms:
    leafset.setdefault(find(z), frozenset())

changed = True
provenance = {find(b): [] for b in chain}
while changed:
    changed = False
    for o, (u, v, prov) in or_children.items():
        ru, rv, ro = find(u), find(v), find(o)
        if ru in leafset and rv in leafset:
            new = leafset[ru] | leafset[rv]
            if ro not in leafset or len(new) > len(leafset[ro]):
                leafset[ro] = new
                provenance[ro] = provenance.get(ru, []) + provenance.get(rv, []) + prov
                changed = True

sizes = Counter(len(v) for v in leafset.values())
print("OR-wire leaf-count histogram:", dict(sorted(sizes.items())))
mx = max(sizes)
roots = [k for k, v in leafset.items() if len(v) == mx]
print(f"largest OR node covers {mx} selectors; {len(roots)} such classes")

# which wires (any alias) are pinned to 1 ?
pinned_class = {}
for w, aid in pin1.items():
    pinned_class.setdefault(find(w), []).append((w, aid))

hits = []
for cls, S in leafset.items():
    if cls in pinned_class and len(S) > 0:
        hits.append((cls, sorted(S), pinned_class[cls]))
hits.sort(key=lambda h: -len(h[1]))
print(f"\nOR wires pinned to 1: {len(hits)}")
for cls, S, pins in hits[:10]:
    eqs = sorted({e for w, aid in pins for e in atom_eqs[aid]})
    print(f"  class rep x_{cls}: OR of {len(S)} selectors = 1   "
          f"pins {[(f'x_{w}', f'atom#{a}') for w, a in pins]} in eqs {eqs[:10]}")

out = {'n_or_gadgets': len(ors),
       'max_or_fanin': mx,
       'or_constraints': [
           {'wire_class': cls, 'n_selectors': len(S), 'selectors': S,
            'pins': [{'wire': w, 'atom': a,
                      'equations': sorted(set(atom_eqs[a]))} for w, a in pins],
            'gadget_atoms': sorted(set(provenance.get(cls, [])))}
           for cls, S, pins in hits]}
json.dump(out, open(os.path.join(HERE, 'ortree.json'), 'w'))
print("wrote ortree.json")
