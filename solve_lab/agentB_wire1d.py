#!/usr/bin/env python3
"""Can the wire SPLIT? For each unpacking eq, find 'lone' wire members (net nonzero linear coeff
after pairing) -> must keep pinned value. For core eqs, find which wire members they use (want!=p).
Then check whether the identity-graph ties the two sets together (forcing a single value)."""
import json, pickle
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
env.forced = {v: (s % p) for v, s in wire.items()}
env.set_from_solution(best)

UNPACK = [8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666]
CORE = [2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]

def wire_terms(i):
    """net linear-wire coeff, and wire members in products, for equation root i."""
    lin = defaultdict(int); prod = set()
    for m, c in env.root_poly[i].items():
        wm = [v for v in m if v in wireset]
        if len(m) == 1 and len(wm) == 1:
            lin[wm[0]] += c
        elif wm:
            prod |= set(wm)
    lone = {v: c for v, c in lin.items() if c % p != 0}
    return lone, prod

print("[d] LONE wire members (net nonzero linear coeff) per unpacking eq:")
lone_all = set()
for i in UNPACK:
    lone, prod = wire_terms(i)
    lone_all |= set(lone)
    print(f"   eq {i}: lone={ {v:c for v,c in lone.items()} }, wire-in-products={sorted(prod)}")
print(f"[d] union of lone (checked) wire members: {sorted(lone_all)}")

print("\n[d] wire members used by CORE eqs:")
core_wire_lin = set(); core_wire_prod = set()
for i in CORE:
    lone, prod = wire_terms(i)
    core_wire_lin |= set(lone); core_wire_prod |= prod
print(f"   core lone-linear wire: {sorted(core_wire_lin)}")
print(f"   core wire-in-products: {sorted(core_wire_prod)}")
print(f"   overlap lone(unpack) ∩ core-wire: {sorted(lone_all & (core_wire_lin|core_wire_prod))}")

# does any core eq need a specific wire member != p? Which wire members multiply handles in core?
# M1 = L1 + x_5101*x_30317 -> x_5101 in product with handle x_30317
handles = {30317, 2936, 5146}
print("\n[d] wire members multiplying the quotient handles in core products:")
for i in CORE:
    for m, c in env.root_poly[i].items():
        wm = [v for v in m if v in wireset]
        hm = [v for v in m if v in handles]
        if wm and hm:
            print(f"   eq {i}: term {c} * wire{wm} * handle{hm}")
            break

# identity graph: are the lone (checked) members in the SAME union-find class as core-product members?
# they are ALL in the 220 class by construction. Question: can we split? Check how many DISTINCT
# lone members vs core members, and whether keeping lone at p / rest at 1 keeps unpacking satisfied.
env.forced = {v: (s % p) for v, s in wire.items()}
for v in lone_all: env.forced[v] = p % p  # keep lone at p, rest at 1
env.set_from_solution(best)
res = env.all_root_residuals()
unpack_fail = [i for i in UNPACK if i in res]
core_fail = [i for i in CORE if i in res]
other_fail = [i for i in res if i not in UNPACK and i not in CORE]
print(f"\n[d] SPLIT wire (lone={len(lone_all)} at p, rest at 1): total fail={len(res)}")
print(f"   unpacking still failing: {len(unpack_fail)} {unpack_fail}")
print(f"   core failing: {len(core_fail)}")
print(f"   OTHER (newly broken by split): {len(other_fail)} {sorted(other_fail)[:20]}")
