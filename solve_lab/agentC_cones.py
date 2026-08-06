#!/usr/bin/env python3
"""Compute free-input cones of the 5 deep vars + S,T; report sizes/overlaps and which free inputs
move each deep var (linear response mod p). Uses partial-forward for speed."""
import json, sys
from agentC_common import (p, gates, order, definer, gcode, forward, partial_forward, downstream_ks,
                           val, freeinp, anc, ns, lines, eqcode, eqvars, load_best, CORE, posof,
                           NVARS, pinned)
from collections import defaultdict

gate_defs = {t: (rhs, vids) for t, rhs, vids in gates}
def freecone(root):
    seen=set(); leaves=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x][1]: st.append(u)
        elif x in freeinp: leaves.add(x)
    return leaves

best = load_best()
DEEP = [3558, 29322, 33469, 27713, 1326]
NAMES = {3558:'x_3558',29322:'x_29322',33469:'x_33469',27713:'x_27713',1326:'x_1326',35389:'S',6671:'T'}
cones = {d: freecone(d) for d in DEEP+[35389,6671]}
for d in DEEP+[35389,6671]:
    print(f"{NAMES[d]} cone: {len(cones[d])} free inputs")
union = set().union(*[cones[d] for d in DEEP])
print(f"union of 5 deep cones: {len(union)} free inputs")
# how many are 'determined' (nonzero in best) vs zero
det = set(v for v in union if best.get(v,0)!=0)
print(f"  of which nonzero-in-best (determined): {len(det)}; zero: {len(union)-len(det)}")

# baseline residues
forward()
base = {d: val[d]%p for d in DEEP+[35389,6671]}
# For speed, precompute downstream ks per handle only within union (bounded)
# measure linear response of each deep var to each handle in union
resp = defaultdict(dict)  # handle -> {deepvar: delta mod p}
movers_of = defaultdict(set)  # deepvar -> set of handles that move it
handles = sorted(union)
print(f"measuring linear response of {len(handles)} handles...", flush=True)
import time
t0=time.time()
for hi,h in enumerate(handles):
    ks = downstream_ks(h)
    o = val[h]; val[h] = o+1
    partial_forward(ks)
    for d in DEEP+[35389,6671]:
        dd = (val[d]%p - base[d])%p
        if dd:
            resp[h][d]=dd
            if d in DEEP: movers_of[d].add(h)
    val[h]=o
    partial_forward(ks)
    if hi%200==0: print(f"  {hi}/{len(handles)} ({time.time()-t0:.0f}s)", flush=True)
print(f"done in {time.time()-t0:.0f}s")
for d in DEEP:
    print(f"movers of {NAMES[d]}: {len(movers_of[d])}")
# handles that move S or T
movers_S = set(h for h in handles if 35389 in resp[h])
movers_T = set(h for h in handles if 6671 in resp[h])
print(f"movers of S: {len(movers_S)}; movers of T: {len(movers_T)}")
# save response data
out = {'handles': handles, 'base': {str(d):base[d] for d in base},
       'resp': {str(h): {str(d):resp[h][d] for d in resp[h]} for h in handles if resp[h]}}
json.dump(out, open('agentC_resp.json','w'))
print("saved agentC_resp.json")
