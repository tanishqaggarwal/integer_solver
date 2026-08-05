#!/usr/bin/env python3
"""Extract and analyze the bounded conflict cluster.

Only ~81 vars are entangled with the 4 violated atoms. Isolate that cluster with
the rest of the circuit fixed to best, find its internal gate orientation and its
true free inputs, and report their nature (bits? residue values?) so we can solve
the small subsystem directly."""
import json, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
NV_PRIM = 4

atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
best = json.load(open('best/best_partial_39013.json'))
ov = json.load(open('cand_override.json'))
bval = [0]*NVARS
for k, x in best.items(): bval[int(k[2:])] = x

changed = set(int(k[2:]) for k in ov if best.get(k, 0) != ov[k])
viol_atoms = [8461, 8463, 20846, 26775, 30862, 31657, 34408, 40873, 43244, 44603, 45293]
cluster = set(changed)
for a in viol_atoms: cluster |= avars[a]
# grow one hop through primitive gates that share >=2 cluster vars (entangled gates)
for _ in range(2):
    add = set()
    for a in range(len(atoms)):
        if len(avars[a]) <= NV_PRIM and len(avars[a] & cluster) >= 2:
            add |= avars[a]
    if add <= cluster: break
    cluster |= add
print(f"cluster size: {len(cluster)}")

# boolean vars (bits): appear in an atom x*(x-1) i.e. {(v,v):1,(v,):-1} or {(v,):..}
# detect via atoms of form v^2 - v
bits = set()
for a in range(len(atoms)):
    poly = atoms[a]
    if len(poly) == 2:
        ms = sorted(poly.keys(), key=len)
        if len(ms[0]) == 1 and len(ms[1]) == 2:
            v = ms[0][0]
            if ms[1] == (v, v) and poly[ms[0]] == -poly[ms[1]]:
                bits.add(v)
cluster_bits = cluster & bits
print(f"cluster bits (boolean): {len(cluster_bits)}: {sorted(cluster_bits)}")

# internal matching (product-priority) among cluster, boundary fixed
cand = defaultdict(list); has_prod = {}; outv = {}
for a in range(len(atoms)):
    if len(avars[a]) > NV_PRIM: continue
    if not (avars[a] & cluster): continue
    bad = set(); lin = set(); hp = False
    for m in atoms[a]:
        if len(m) == 1:
            if m[0] in cluster: lin.add(m[0])
        else:
            hp = True
            for x in m:
                if x in cluster: bad.add(x)
    outv[a] = lin - bad; has_prod[a] = hp
    for v in (lin - bad): cand[v].append(a)
match = {}; used = {}
for a in outv:
    if has_prod[a] and len(outv[a]) == 1:
        v = next(iter(outv[a]))
        if v not in match and a not in used: match[v] = a; used[a] = v
sys.setrecursionlimit(100000)
def aug(v, seen):
    for a in cand[v]:
        if a in seen: continue
        seen.add(a)
        if a not in used:
            match[v] = a; used[a] = v; return True
        if has_prod[a]: continue
        if aug(used[a], seen):
            match[v] = a; used[a] = v; return True
    return False
for v in list(cluster):
    if v not in match and cand[v]: aug(v, set())
free = sorted(v for v in cluster if v not in match)
print(f"cluster internal: matched {len(match)}, free inputs {len(free)}")
print(f"free inputs: {free}")
print(f"free that are bits: {sorted(set(free) & bits)}")
# for each free input, is best value 0/1 or big?
for v in free:
    bv = bval[v]
    tag = 'BIT' if v in bits else ('0/1' if bv in (0,1) else 'BIG' if abs(bv) > 10**50 else str(bv))
    print(f"  x_{v}: best={tag}" + ('' if tag not in ('BIG',) else f' ({len(str(abs(bv)))} digits)'))
