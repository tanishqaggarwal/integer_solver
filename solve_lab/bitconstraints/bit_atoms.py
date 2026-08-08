#!/usr/bin/env python3
"""Pass 10: EXHAUSTIVE classification of every atom that mentions a selector bit.

If every such atom falls into a known template, then the only way the selectors
influence the rest of the system is through those templates -- which pins down
exactly what the system can say about the bits alone.
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

touch = [i for i, a in enumerate(atoms) if support(a) & BITS]
print(f"atoms mentioning >=1 selector bit: {len(touch)}")

cls = defaultdict(list)
for i in touch:
    a = atoms[i]
    d = dict(a)
    s = support(a)
    sb = s & BITS
    so = s - BITS
    dg = degree(a)
    tag = None
    if len(s) == 1 and dg == 2 and len(a) == 2:
        (v,) = s
        if d.get((v, v)) == -d.get((v,)) != 0:
            tag = 'booleanity b - b^2'
    if tag is None and len(s) == 2 and dg == 1 and len(a) == 2 and () not in d:
        tag = 'copy  b - w'
    if tag is None and len(s) == 2 and dg == 1 and len(a) == 3 and d.get(()) == 1:
        u, v = sorted(s)
        if d[(u,)] == -1 and d[(v,)] == -1:
            tag = 'NOT   1 - b - w'
    if tag is None and dg == 2 and len(sb) >= 1 and all(len(m) != 1 or m[0] not in BITS for m in d):
        pass
    if tag is None:
        # gated pin:  b*(w - C)  ==  b*w - C*b        (2 terms, deg 2, 2 vars)
        if len(a) == 2 and dg == 2 and len(s) == 2 and len(sb) == 1:
            (b,) = sb
            (w,) = so
            if d.get(tuple(sorted((b, w)))) is not None and d.get((b,)) is not None:
                tag = 'pin-gate  b*(w - C)'
    if tag is None:
        # complement gate:  (1-b)*w  ==  w - b*w   (2 terms deg2, 2 vars)
        if len(a) == 2 and dg == 2 and len(s) == 2 and len(sb) == 1:
            (b,) = sb
            (w,) = so
            if d.get((w,)) is not None and d.get(tuple(sorted((b, w)))) is not None:
                tag = 'complement-gate (1-b)*w'
    if tag is None:
        # sums of booleanity terms only?
        if all(len(m) in (1, 2) and (len(m) == 1 or m[0] == m[1]) for m in d) and () not in d:
            byv = defaultdict(lambda: [0, 0])
            for m, c in d.items():
                byv[m[0]][len(m) - 1] += c
            if all(x[0] == -x[1] for x in byv.values()):
                tag = 'sum of booleanity terms'
    if tag is None:
        tag = f'OTHER nb={len(sb)} nother={len(so)} deg={dg} nterms={len(a)}'
    cls[tag].append(i)

print("\nclassification:")
for k in sorted(cls, key=lambda k: -len(cls[k])):
    print(f"  {len(cls[k]):6d}  {k}")

others = [i for k, v in cls.items() if k.startswith('OTHER') for i in v]
print(f"\nunclassified atoms mentioning a bit: {len(others)}")
for i in others[:60]:
    a = atoms[i]
    print(f"  atom#{i} bits={sorted(support(a)&BITS)} others={len(support(a)-BITS)} "
          f"deg={degree(a)} nterms={len(a)} eqs={sorted(set(atom_eqs[i]))[:5]}")
    print(f"      {fmt(a)[:260]}")

# --- per-bit inventory ---
per = defaultdict(Counter)
for k, v in cls.items():
    base = k.split()[0]
    for i in v:
        for b in support(atoms[i]) & BITS:
            per[b][base] += 1
sig = Counter(tuple(sorted(per[b].items())) for b in BITS)
print("\nper-selector-bit atom signature (identical for all bits?):")
for k, v in sig.most_common():
    print(f"   x{v} bits: {k}")

json.dump({k: v for k, v in cls.items()},
          open(os.path.join(HERE, 'bit_atom_classes.json'), 'w'))
print("wrote bit_atom_classes.json")
