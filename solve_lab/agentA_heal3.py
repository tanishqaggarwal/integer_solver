#!/usr/bin/env python3
"""Elegant heal: align x_19083 = x_24908 via free input x_21589 in its definition,
so x_16742 = x_19083 = x_24908 (all checks hold + x_3558=0). Same for x_12186 side."""
import json
from agentA_harness import (p, load_solution, forward, gates, definer, order, freeinp,
                            backward_cone, eqcode, eqvars, NEQ, NVARS)

base = load_solution('best/new_instance_partial_39013.json'); forward(base)
def count(v):
    ns = {'__builtins__': {}, 'v': v}
    return set(i for i in range(NEQ) if eval(eqcode[i], ns) != 0)
F0 = count(base)

gdef = {t: gates[definer[t]][1] for t in order}
gvids = {t: gates[definer[t]][2] for t in order}

# check x_21589 controls x_19083 and is independent of x_24908
_, c24908 = backward_cone(24908)
_, c19083 = backward_cone(19083)
print(f"x_21589 in x_24908 cone? {21589 in c24908}  (must be False)")
print(f"x_21589 in x_19083 cone? {21589 in c19083}")
print(f"x_19083 def: {gdef[19083]}   x_23758 def: {gdef.get(23758)}")
print(f"x_38170 def: {gdef.get(38170)}   value mod p = {base[38170]%p}  exact={base[38170]}")
print(f"x_6361={base[6361]}  (mod p={base[6361]%p})")
print(f"x_24908 - x_6361 = {base[24908]-base[6361]}")
print(f"x_38170 divides (x_24908 - x_6361)? {(base[24908]-base[6361]) % base[38170] == 0 if base[38170]!=0 else 'x_38170=0'}")

# x_19083 = x_6361 + x_38170*x_21589 ; solve x_21589 for x_19083 = x_24908
if base[38170] != 0 and (base[24908]-base[6361]) % base[38170] == 0:
    x21589_new = (base[24908]-base[6361]) // base[38170]
    print(f"exact x_21589 = {x21589_new}")
else:
    x21589_new = None
    print("no exact integer x_21589 for x_19083=x_24908")

# Now the x_12186 side: what pins x_12186? find its difference partner gate
import re
lines = [L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
dp = re.compile(r'\(x_12186\)\s*-\s*\(x_(\d+)\)'); dp2 = re.compile(r'\(x_(\d+)\)\s*-\s*\(x_12186\)')
from collections import Counter
pa, pb = Counter(), Counter()
for i in range(NEQ):
    if 12186 in eqvars[i]:
        for m in dp.finditer(lines[i]): pa[int(m.group(1))]+=1
        for m in dp2.finditer(lines[i]): pb[int(m.group(1))]+=1
print(f"\nx_12186 partners (x_12186 - x_Y): {dict(pa)}")
print(f"x_12186 partners (x_Y - x_12186): {dict(pb)}")
for y in list(pa)+list(pb):
    print(f"  x_{y}: free={y in freeinp} def={gdef.get(y,'FREE')[:50]} base==x12186? {base[y]==base[12186]}")

# Build A2 + alignment and test
def build(align19083):
    v = base[:]
    if align19083 and x21589_new is not None:
        v[21589] = x21589_new
    v[16742] = base[24908]           # x_3558 = 0 (x_16742 = x_24908 = new x_19083)
    v[12186] = base[14853]           # x_29322 = 0 (x_12186 = x_14853 = x_1308)
    forward(v)
    if v[11150] % p == 0: v[30317] = -(v[11150])//p
    if (537773*v[37758]) % p == 0: v[2936] = (537773*v[37758])//p
    if v[25739] % (6672769*p) == 0: v[5146] = v[25739]//(6672769*p)
    return v

for al in [False, True]:
    v = build(al); F = count(v)
    print(f"\nalign_x19083={al}: {NEQ-len(F)}/{NEQ} ({len(F)} fail)")
    print(f"  x_19083==x_24908? {v[19083]==v[24908]}  x_16742==x_19083? {v[16742]==v[19083]}")
    print(f"  S,T mod p = {v[35389]%p==0},{v[6671]%p==0}")
    broke = sorted(F-F0); print(f"  broke {len(broke)}: {broke}")
    if NEQ-len(F) > 39016:
        json.dump({f"x_{i}":v[i] for i in range(NVARS) if v[i]!=0}, open(f'best_agentA_{NEQ-len(F)}.json','w'))
        print(f"  SAVED best_agentA_{NEQ-len(F)}.json")
