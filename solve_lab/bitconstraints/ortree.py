#!/usr/bin/env python3
"""Pass 7: reconstruct the boolean layer over the selector bits.

Wire vocabulary recovered from atom shapes:
    pin      : 1 - x            ->  x = 1
    zero     : x                ->  x = 0
    copy     : x - y            ->  x = y
    not      : 1 - x - y        ->  y = 1 - x
    and      : x*y - z          ->  z = x AND y
Every wire reachable from the 256 selectors through these gates is a boolean
function of the bits of the restricted form  OR(S)  or  NOR(S) = AND_{i in S}(1-b_i).
We compute that symbolic form for every reachable wire, then look for wires that
are pinned to a constant -- those are the pure-bit cardinality constraints.
"""
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

pin_one, pin_zero, copies, nots, ands, others2 = [], [], [], [], [], []
for i, a in enumerate(atoms):
    d = dict(a)
    s = support(a)
    dg = degree(a)
    if len(a) == 1 and len(s) == 1 and dg == 1 and d[(next(iter(s)),)] == 1:
        pin_zero.append((i, next(iter(s))))
    elif len(a) == 2 and len(s) == 1 and dg == 1 and () in d:
        v = next(iter(s))
        c0, c1 = d[()], d[(v,)]
        if c0 == 1 and c1 == -1:
            pin_one.append((i, v))
        else:
            others2.append((i, v, c0, c1))
    elif len(a) == 2 and len(s) == 2 and dg == 1:
        u, v = sorted(s)
        if d[(u,)] == 1 and d[(v,)] == -1:
            copies.append((i, u, v))
    elif len(a) == 3 and len(s) == 2 and dg == 1 and () in d:
        u, v = sorted(s)
        if d[()] == 1 and d[(u,)] == -1 and d[(v,)] == -1:
            nots.append((i, u, v))          # u + v = 1
    elif len(a) == 2 and len(s) == 3 and dg == 2:
        q = [m for m in d if len(m) == 2]
        l = [m for m in d if len(m) == 1]
        if len(q) == 1 and len(l) == 1 and d[q[0]] == -d[l[0]] and q[0][0] != q[0][1]:
            ands.append((i, q[0][0], q[0][1], l[0][0]))

print(f"pins  x = 1        : {len(pin_one)}")
print(f"pins  x = 0        : {len(pin_zero)}")
print(f"copies x = y       : {len(copies)}")
print(f"NOT   x + y = 1    : {len(nots)}")
print(f"AND   x*y = z      : {len(ands)}")
print(f"other 1-var linear : {len(others2)}  "
      f"coef-pairs {Counter((c0,c1) for _,_,c0,c1 in others2).most_common(6)}")

# ---------------- symbolic boolean evaluation -------------------------
# form: ('OR', frozenset) or ('NOR', frozenset).  bit b  == OR({b}).
def neg(f):
    return ('NOR' if f[0] == 'OR' else 'OR', f[1])


form = {b: ('OR', frozenset([b])) for b in chain}
prov = {b: [] for b in chain}

adj_copy = defaultdict(list)
for i, u, v in copies:
    adj_copy[u].append((v, i))
    adj_copy[v].append((u, i))
adj_not = defaultdict(list)
for i, u, v in nots:
    adj_not[u].append((v, i))
    adj_not[v].append((u, i))
and_by_in = defaultdict(list)
for i, u, v, w in ands:
    and_by_in[u].append((i, u, v, w))
    and_by_in[v].append((i, u, v, w))

frontier = list(chain)
while frontier:
    nxt = []
    for x in frontier:
        for y, i in adj_copy[x]:
            if y not in form:
                form[y] = form[x]; prov[y] = prov[x] + [i]; nxt.append(y)
        for y, i in adj_not[x]:
            if y not in form:
                form[y] = neg(form[x]); prov[y] = prov[x] + [i]; nxt.append(y)
        for i, u, v, w in and_by_in[x]:
            if w in form or u not in form or v not in form:
                continue
            fu, fv = form[u], form[v]
            if fu[0] == 'NOR' and fv[0] == 'NOR':
                form[w] = ('NOR', fu[1] | fv[1])
            elif fu[0] == 'OR' and fv[0] == 'OR' and fu[1] == fv[1]:
                form[w] = fu                      # idempotent
            else:
                continue
            prov[w] = prov[u] + prov[v] + [i]
            nxt.append(w)
    frontier = nxt

print(f"\nwires with a recovered OR/NOR form: {len(form)}")
hist = Counter((f[0], len(f[1])) for f in form.values())
print("form histogram (kind, |S|) -> count:")
for k in sorted(hist, key=lambda k: (k[0], k[1])):
    print(f"   {k} -> {hist[k]}")

full = [v for v, f in form.items() if len(f[1]) == 256]
print(f"\nwires depending on ALL 256 bits: {len(full)}  "
      f"kinds {Counter(form[v][0] for v in full)}")

pin1 = {v: i for i, v in pin_one}
pin0 = {v: i for i, v in pin_zero}
constraints = []
for v, f in form.items():
    if v in pin1:
        constraints.append((v, f[0], sorted(f[1]), 1, pin1[v]))
    if v in pin0:
        constraints.append((v, f[0], sorted(f[1]), 0, pin0[v]))
print(f"\npinned OR/NOR wires (=> cardinality constraints on the bits): "
      f"{len(constraints)}")
seen = set()
out = []
for v, kind, S, val, aid in constraints:
    key = (kind, tuple(S), val)
    eqs = sorted(set(atom_eqs[aid]))
    txt = (f"{kind}({len(S)} bits) = {val}")
    print(f"  x_{v}: {txt}   via atom#{aid}, equations {eqs[:12]}"
          f"{' ...' if len(eqs)>12 else ''}")
    if key not in seen:
        seen.add(key)
        out.append({'wire': v, 'kind': kind, 'bits': S, 'value': val,
                    'atom': aid, 'equations': eqs,
                    'gate_chain_atoms': prov[v]})
json.dump(out, open(os.path.join(HERE, 'cardinality.json'), 'w'))
json.dump({'pin_one': pin_one, 'pin_zero': pin_zero,
           'n_copies': len(copies), 'n_nots': len(nots), 'n_ands': len(ands)},
          open(os.path.join(HERE, 'boolean_layer.json'), 'w'))
print("wrote cardinality.json, boolean_layer.json")
