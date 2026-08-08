#!/usr/bin/env python3
"""Pass 6: global atom shape census + explicit hunt for the OR-tree over the
selector bits and for unit pins (x = 1)."""
import os, json
from collections import Counter, defaultdict
from scan import load, support, degree, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain, BITS = load()
atoms = D['atoms']
atom_eqs = defaultdict(list)
for i, terms in enumerate(D['eq_terms']):
    for c, aid in terms:
        atom_eqs[aid].append(i)


def sig(poly):
    """shape signature: sorted (monomial-arity, coef) pairs"""
    return tuple(sorted((len(m), c) for m, c in poly))


# --- unit pins  x - 1 ---
pins1 = [i for i, a in enumerate(atoms) if len(a) == 2 and
         dict(a).get(()) is not None and len(support(a)) == 1 and degree(a) == 1]
print(f"linear 2-term atoms with a constant (pins/NOTs): {len(pins1)}")
pin_one, pin_other = [], []
for i in pins1:
    d = dict(atoms[i])
    (v,) = support(atoms[i])
    if d.get(()) == -1 and d.get((v,)) == 1:
        pin_one.append((i, v))
    else:
        pin_other.append((i, v, d))
print(f"  exact 'x - 1' pins: {len(pin_one)}")
print(f"  other c0 + c1*x   : {len(pin_other)}  sample {pin_other[:5]}")

# --- OR gates:  a + b - a*b - c   (c = a OR b) ---
ors = []
for i, a in enumerate(atoms):
    if len(a) != 4 or degree(a) != 2:
        continue
    s = support(a)
    if len(s) != 3:
        continue
    d = dict(a)
    quad = [m for m in d if len(m) == 2]
    lin = [m for m in d if len(m) == 1]
    if len(quad) != 1 or len(lin) != 3:
        continue
    (u, v) = quad[0]
    if u == v:
        continue
    w = (s - {u, v})
    if len(w) != 1:
        continue
    w = w.pop()
    cq, cu, cv, cw = d[quad[0]], d[(u,)], d[(v,)], d[(w,)]
    # a+b-ab-c  (up to overall sign)
    if cu == cv == -cq and cw == -cu:
        ors.append((i, u, v, w))
print(f"\nOR/AND-complement gates  (a + b - a*b - c): {len(ors)}")

# --- AND gates: a*b - c ---
ands = []
for i, a in enumerate(atoms):
    if len(a) != 2 or degree(a) != 2:
        continue
    s = support(a)
    if len(s) != 3:
        continue
    d = dict(a)
    q = [m for m in d if len(m) == 2]
    l = [m for m in d if len(m) == 1]
    if len(q) == 1 and len(l) == 1 and d[q[0]] == -d[l[0]]:
        ands.append((i, q[0][0], q[0][1], l[0][0]))
print(f"AND gates (a*b - c): {len(ands)}")

# --- copies / NOTs touching bits ---
copies = [(i, sorted(support(a))) for i, a in enumerate(atoms)
          if len(a) == 2 and degree(a) == 1 and len(support(a)) == 2
          and sorted(c for m, c in a) == [-1, 1]]
nots = [(i, sorted(support(a))) for i, a in enumerate(atoms)
        if len(a) == 3 and degree(a) == 1 and len(support(a)) == 2
        and dict(a).get(()) is not None]
print(f"copy atoms (x - y): {len(copies)}")
print(f"NOT-ish atoms (c + a*x + b*y): {len(nots)}")

# --- build the boolean gate graph and find the OR-tree over the bits ---
# node -> expression over bits, as a frozenset of bits it ORs together
orset = {b: frozenset([b]) for b in BITS}
gate = {}
for i, u, v, w in ors:
    gate.setdefault(w, []).append(('or', u, v, i))
for i, u, v, w in ands:
    gate.setdefault(w, []).append(('and', u, v, i))
# copies propagate
copy_of = defaultdict(list)
for i, (a, b) in copies:
    copy_of[a].append((b, i))
    copy_of[b].append((a, i))

changed = True
prov = {b: [] for b in BITS}
while changed:
    changed = False
    for w, gs in gate.items():
        for kind, u, v, i in gs:
            if kind == 'or' and u in orset and v in orset:
                new = orset[u] | orset[v]
                if orset.get(w) != new and (w not in orset or len(new) > len(orset[w])):
                    orset[w] = new
                    prov[w] = prov.get(u, []) + prov.get(v, []) + [i]
                    changed = True
    for a, lst in list(copy_of.items()):
        if a in orset:
            for b, i in lst:
                if b not in orset or len(orset[a]) > len(orset.get(b, ())):
                    orset[b] = orset[a]
                    prov[b] = prov.get(a, []) + [i]
                    changed = True

best = max(orset.items(), key=lambda kv: len(kv[1]))
print(f"\nlargest OR-of-bits wire: x_{best[0]} = OR of {len(best[1])} selector bits")
sizes = Counter(len(v) for v in orset.values())
print("OR-fanin size histogram:", dict(sorted(sizes.items())[:12]),
      "... max", max(sizes))
bigs = sorted([k for k, v in orset.items() if len(v) == len(best[1])])
print(f"wires that OR all {len(best[1])} bits: {len(bigs)} -> {bigs[:20]}")

pinned_ones = {v for i, v in pin_one}
hit = [v for v in bigs if v in pinned_ones]
print(f"...of which are pinned to 1 by an 'x - 1' atom: {len(hit)} -> {hit[:20]}")
for v in hit[:10]:
    pa = [i for i, w in pin_one if w == v]
    print(f"   x_{v} = 1 via atom(s) {pa}, in equations "
          f"{sorted({e for i in pa for e in atom_eqs[i]})[:10]}")

json.dump({
    'or_gates': ors, 'and_gates': ands,
    'unit_pins': pin_one,
    'or_all_wires': bigs,
    'or_all_pinned_to_one': hit,
    'or_all_size': len(best[1]),
}, open(os.path.join(HERE, 'gates.json'), 'w'))
print("wrote gates.json")
