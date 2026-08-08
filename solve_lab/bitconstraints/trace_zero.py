#!/usr/bin/env python3
"""Pass 9: instrumented propagation -- record which atom determined each wire,
then back-trace the all-zeros contradiction to the responsible gate chain."""
import os, json, sys
from collections import defaultdict
from closure import load, reduce_poly, atom_vars, NVARS
from scan import fmt

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain = load()
BITS = set(chain)
atoms = D['atoms']
atom_eqs = defaultdict(list)
for i, terms in enumerate(D['eq_terms']):
    for c, aid in terms:
        atom_eqs[aid].append(i)


def propagate_traced(val):
    watch = defaultdict(list)
    for i, a in enumerate(atoms):
        for v in atom_vars(a):
            watch[v].append(i)
    src = {}                        # var -> (atom_id, [input vars])
    queue = list(range(len(atoms)))
    inq = [True]*len(atoms)
    contra = []
    while queue:
        i = queue.pop(); inq[i] = False
        red = reduce_poly(atoms[i], val)
        uv = set()
        for m in red:
            uv.update(m)
        if not uv:
            if red:
                contra.append((i, red.get((), 0)))
            continue
        if len(uv) != 1:
            continue
        u = next(iter(uv))
        c = [0, 0, 0, 0]
        bad = False
        for m, k in red.items():
            if len(m) > 3:
                bad = True; break
            c[len(m)] += k
        if bad:
            continue
        nv = None
        if c[2] == 0 and c[3] == 0:
            if c[1] == 0:
                if c[0] != 0: contra.append((i, c[0]))
                continue
            if c[0] % c[1]: contra.append((i, 'nondiv')); continue
            nv = -c[0]//c[1]
        elif c[3] == 0:
            disc = c[1]*c[1]-4*c[2]*c[0]
            if disc < 0: contra.append((i, 'negdisc')); continue
            r = int(disc**0.5)
            while r*r > disc: r -= 1
            while (r+1)*(r+1) <= disc: r += 1
            if r*r != disc: contra.append((i, 'nonsq')); continue
            roots = {(-c[1]+s)//(2*c[2]) for s in (r, -r)
                     if (-c[1]+s) % (2*c[2]) == 0}
            if len(roots) != 1: continue
            nv = roots.pop()
        else:
            continue
        val[u] = nv
        src[u] = (i, sorted(atom_vars(atoms[i]) - {u}))
        for j in watch[u]:
            if not inq[j]:
                inq[j] = True; queue.append(j)
    return src, contra


val = [None]*NVARS
for b in chain:
    val[b] = 0
src, contra = propagate_traced(val)
print(f"determined {sum(1 for v in val if v is not None)}; "
      f"contradictions {len(contra)}: {contra}")

for aid, resid in contra:
    print(f"\n=== contradicting atom#{aid} (residual {resid}) ===")
    print("   ", fmt(atoms[aid])[:300])
    print("    appears in equations:", sorted(set(atom_eqs[aid]))[:12])
    # back-trace each variable of the atom
    seen = set()
    order = []

    def walk(v, depth):
        if v in seen or depth > 40:
            return
        seen.add(v)
        if v in BITS:
            order.append((v, 'SELECTOR', val[v], None)); return
        if v not in src:
            order.append((v, 'FREE/UNDET', val[v], None)); return
        a, ins = src[v]
        order.append((v, 'atom#%d' % a, val[v], ins))
        for u in ins:
            walk(u, depth+1)

    for v in sorted(atom_vars(atoms[aid])):
        walk(v, 0)
    print(f"    back-cone size {len(order)}; "
          f"selectors in cone: {sum(1 for o in order if o[1]=='SELECTOR')}; "
          f"undetermined in cone: {sum(1 for o in order if o[1]=='FREE/UNDET')}")
    for v, how, value, ins in order[:40]:
        print(f"      x_{v} = {str(value)[:40]}  <- {how}  from {ins if ins is None or len(ins)<8 else str(ins[:8])+'...'}")
