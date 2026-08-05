#!/usr/bin/env python3
"""Extract the wiring null space and characterize the obstruction precisely.
1. rank(J_sat) and its null space (over active columns).
2. For each of the 20 core roots: directional derivative along each null vector (immovable=0).
3. Rank of the 20 core-root gradients MODULO the wiring row space (# independent core conditions),
   and whether the residuals are consistent.
4. Finite move along the null vector(s): does it change stuck core roots at 2nd order?"""
import json, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best); env.tangent_linear()
res0 = env.all_root_residuals(); core = sorted(res0)

# active columns
active = set()
sat_rows = []; core_rows = []
for i in range(len(env.root_poly)):
    g = env.root_grad(i)
    if not g: continue
    for c in g: active.add(c)
    (core_rows if i in res0 else sat_rows).append((i, g))
active = sorted(active)
aidx = {c: j for j, c in enumerate(active)}
print(f"[null] active cols={len(active)}, sat rows={len(sat_rows)}, core rows={len(core_rows)}")

# RREF of J_sat over active columns; record pivot columns and reduced pivot rows
def rref(rows):
    pivots = {}  # col -> rowdict (normalized, pivot=1), only over active cols
    coldeg = defaultdict(int)
    for _, rd in rows:
        for c in rd: coldeg[c] += 1
    order = sorted(range(len(rows)), key=lambda k: len(rows[k][1]))
    for k in order:
        rd = dict(rows[k][1])
        while True:
            pcol = None
            for c in rd:
                if c in pivots: pcol = c; break
            if pcol is None: break
            f = rd[pcol]; prow = pivots[pcol]
            for c, v in prow.items():
                nv = (rd.get(c, 0) - f*v) % p
                if nv: rd[c] = nv
                elif c in rd: del rd[c]
        if not rd: continue
        pc = min(rd, key=lambda c: coldeg.get(c, 0))
        inv = pow(rd[pc], p-2, p)
        pivots[pc] = {c:(v*inv)%p for c,v in rd.items()}
    return pivots

piv = rref(sat_rows)
rank_sat = len(piv)
pivcols = set(piv)
freecols = [c for c in active if c not in pivcols]
print(f"[null] rank(J_sat)={rank_sat}; nullity over active = {len(active)-rank_sat}; free active cols={freecols}")

# build null vectors: for each free active col f, set f=1, others free=0, back-substitute pivots
def null_vector(freecol):
    n = {freecol: 1}
    # pivot col value = -sum_{c!=pivot} prow[c]*n[c]
    for pc, prow in piv.items():
        s = 0
        for c, v in prow.items():
            if c != pc: s = (s - v * n.get(c, 0)) % p
        n[pc] = s % p
    return {c: v for c, v in n.items() if v}

nulls = [null_vector(f) for f in freecols]
for f, n in zip(freecols, nulls):
    # verify it's in null space of a few sat rows
    ok = True
    for i, g in sat_rows[:200]:
        s = sum(g.get(c, 0) * n.get(c, 0) for c in g) % p
        if s: ok = False; break
    print(f"[null] null vector from free col {f}: support={len(n)}, verified(sat rows)={ok}")

# directional derivative of each core root along each null vector
print("[null] core-root directional derivative along wiring null space:")
movable = 0
for i, g in core_rows:
    ders = []
    for n in nulls:
        d = sum(g.get(c, 0) * n.get(c, 0) for c in g) % p
        ders.append(d)
    nz = any(d != 0 for d in ders)
    if nz: movable += 1
    print(f"  eq {i}: d/dn = {[ '0' if d==0 else 'NONZERO' for d in ders]}  residual={res0[i]!=0}")
print(f"[null] core roots MOVABLE along wiring null space: {movable}/{len(core_rows)}")

# rank of 20 core gradients modulo wiring row space, and consistency
def reduce_against(pivots, rd0, rhs0):
    rd = dict(rd0); rhs = rhs0
    changed = True
    while changed:
        changed = False
        for c in list(rd):
            if c in pivots:
                f = rd[c]; prow = pivots[c]
                for cc, v in prow.items():
                    nv = (rd.get(cc, 0) - f*v) % p
                    if nv: rd[cc] = nv
                    elif cc in rd: del rd[cc]
                changed = True
                break
    return rd, rhs
piv2 = dict(piv); indep=0; incons=0
for i, g in core_rows:
    rd, rhs = reduce_against(piv2, g, (-res0[i]) % p)
    if not rd:
        if rhs % p: incons += 1
    else:
        pc = min(rd); inv = pow(rd[pc], p-2, p)
        piv2[pc] = {c:(v*inv)%p for c,v in rd.items()}; indep += 1
print(f"[null] core modulo wiring: NEW independent LHS directions={indep}, contradictory(0=nonzero)={incons}")

# finite move along the null vector: can 2nd-order change stuck roots?
print("[null] finite move along null vector n (checking 2nd-order effect on core + wiring):")
base = {v: env.valp[v] for v in env.freeset}
n = nulls[0] if nulls else {}
for t in [1, 1000, (p-1)//2, p-1]:
    for v in env.freeset: env.valp[v] = base[v]
    for c, v in n.items():
        env.valp[c] = (env.valp[c] + t * v) % p
    env.forward()
    rv = {i: env.root_val(i) for i in range(len(env.root_poly)) if env.root_val(i)}
    core_change = sum(1 for i in core if rv.get(i, 0) != res0[i])
    wbreak = len(set(rv) - set(core))
    print(f"  t={t if t<10**6 else 'big'}: total fail={len(rv)}, wiring breaks={wbreak}, core changed={core_change}")
for v in env.freeset: env.valp[v] = base[v]
env.forward()
