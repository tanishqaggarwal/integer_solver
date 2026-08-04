#!/usr/bin/env python3
"""Examine best_agentA_39021.json: failing equations (exact Z), wire member values, handle values,
and L2 mod 6672769. Determine what's failing and the wire state."""
import json, pickle
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
sol = {int(k[2:]): int(v) for k, v in json.load(open('best_agentA_39021.json')).items()}
valz = [0]*NVARS
for k, v in sol.items(): valz[k] = v
def root_int(i):
    s = 0
    for m, c in env.root_poly[i].items():
        t = c
        for v in m: t *= valz[v]
        s += t
    return s
fails = [i for i in range(len(env.root_poly)) if root_int(i) != 0]
print(f"[examA] best_agentA_39021: failing={len(fails)}: {fails}")
UNPACK = [8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666]
CORE = [2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]
print(f"[examA] failing that are UNPACK: {sorted(set(fails)&set(UNPACK))}")
print(f"[examA] failing that are CORE:   {sorted(set(fails)&set(CORE))}")
print(f"[examA] failing OTHER: {sorted(set(fails)-set(UNPACK)-set(CORE))}")

# wire values
wv = {v: valz[v] for v in wire}
distinct = set(wv.values())
print(f"[examA] wire member value distinct count: {len(distinct)}; sample:")
for v in [26064, 5101, 38100, 32017, 26789]:
    val = valz[v]
    print(f"   x_{v} = {val}  (=1? {val==1}, =p? {val==p}, mod p={val%p if val%p<10**6 else '...'})")
# how many wire members ==1, ==p, other
c1 = sum(1 for v in wire if valz[v]==1); cp = sum(1 for v in wire if valz[v]==p)
print(f"[examA] wire members ==1: {c1}, ==p: {cp}, other: {len(wire)-c1-cp}")

# handles
for h in [30317, 2936, 5146]:
    print(f"[examA] handle x_{h} = {valz[h]}")

# L2 mod 6672769 : L2 = x_25739 (load 2). print value + mod
for name, v in [("L1=x_11150",11150),("L2=x_25739",25739),("L3=x_37758",37758)]:
    print(f"[examA] {name} = {valz[v]}  ; mod p={valz[v]%p if valz[v]%p<10**7 else '...'} ; mod 6672769={valz[v]%6672769}")
# also x_3558, x_29322
for name, expr in [("x_3558", 3558), ("x_29322", 29322)]:
    print(f"[examA] {name} = {valz[expr]} ; mod p = {valz[expr]%p if valz[expr]%p<10**7 else '...'}")
