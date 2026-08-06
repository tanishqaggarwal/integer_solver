#!/usr/bin/env python3
"""Are the 220 wire members ONE rigid class or splittable? Distinguish HARD ties (equations whose
ROOT is exactly a 2-term wire identity x_a +/- x_b) from SOFT ties (identity only inside larger
combos, healable by handles). Build the hard-tie graph over wire members, find components, and
locate the core-product members vs the lone-unpacking members."""
import json, pickle
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
UNPACK = [8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666]
COREW = {5101, 32017, 26789}

# hard ties: equations whose root poly is exactly {(a,):c, (b,):-c} with a,b wire
hard = defaultdict(set)
nhard = 0
for i, rp in enumerate(env.root_poly):
    if len(rp) != 2: continue
    ms = list(rp.items())
    (m1, c1), (m2, c2) = ms
    if len(m1) == 1 and len(m2) == 1 and m1[0] in wireset and m2[0] in wireset:
        hard[m1[0]].add(m2[0]); hard[m2[0]].add(m1[0]); nhard += 1
print(f"[hard] standalone 2-term wire-identity equations: {nhard}")
# components of hard graph over all 220 wire members
seen = set(); comps = []
for v in wire:
    if v in seen: continue
    stack = [v]; comp = set()
    while stack:
        x = stack.pop()
        if x in seen: continue
        seen.add(x); comp.add(x)
        for u in hard.get(x, ()): stack.append(u)
    comps.append(comp)
comps.sort(key=len, reverse=True)
print(f"[hard] hard-tie components among 220 wire: {len(comps)}; sizes {[len(c) for c in comps][:15]}")
# where are core members and lone unpacking members?
# lone unpacking members = wire members with net-nonzero linear coeff in some unpacking eq
lone = set()
for i in UNPACK:
    lin = defaultdict(int)
    for m, c in env.root_poly[i].items():
        if len(m) == 1 and m[0] in wireset: lin[m[0]] += c
    for v, c in lin.items():
        if c % p: lone.add(v)
for ci, comp in enumerate(comps[:12]):
    cw = comp & COREW; lo = comp & lone
    print(f"  comp {ci} size {len(comp)}: core-members={sorted(cw)}, #lone-unpack={len(lo)}")
print(f"[hard] total lone-unpack members: {len(lone)}; in same comp as core? "
      f"{any((c & COREW) and (c & lone) for c in comps)}")

# also: how many equations reference wire members in NON-identity (soft) ways = ripple potential
soft_eqs = 0
for i, rp in enumerate(env.root_poly):
    if any(v in wireset for m in rp for v in m):
        if not (len(rp) == 2 and all(len(m)==1 and m[0] in wireset for m in rp)):
            soft_eqs += 1
print(f"[hard] equations with wire in non-pure-identity form (ripple-carrying): {soft_eqs}")
