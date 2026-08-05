#!/usr/bin/env python3
"""Examine the exact cones of the small deep vars, and which equations the 97 handles touch."""
import json
from agentC_common import (p, gates, order, definer, gcode, forward, val, freeinp, anc, ns,
                           lines, eqcode, eqvars, load_best, CORE, NVARS, pinned)
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
def show_def(t, depth=0, maxd=3):
    if depth>maxd: return
    if t in gate_defs:
        rhs,vids = gate_defs[t]
        print("  "*depth + f"x_{t} = {rhs}")
        for u in sorted(set(vids)):
            if u in gate_defs and depth<maxd: show_def(u, depth+1, maxd)
            elif u in freeinp: print("  "*(depth+1)+f"x_{u} [FREE, val={val[u]}]")
    else:
        print("  "*depth + f"x_{t} [FREE, val={val[t]}]")

best = load_best()
forward()
for d in [29322, 33469, 27713, 1326]:
    print(f"===== x_{d} (residue {val[d]%p}) =====")
    print("cone free inputs:", sorted(freecone(d)))
    show_def(d, maxd=4)
    print()

# For x_3558: it's x_24908 - x_16742. Show top level
print("===== x_3558 =====")
print("cone size:", len(freecone(3558)))
show_def(3558, maxd=1)
print("x_24908 sub-cone size:", len(freecone(24908)))

# Which currently-satisfied equations do the 97 handles touch?
resp = json.load(open('agentC_resp.json'))
handles = resp['handles']
Hset = set(handles)
# baseline failing set
F0 = set(i for i in range(len(lines)) if eval(eqcode[i], ns)!=0)
touched = set()
for i in range(len(lines)):
    if eqvars[i] & Hset:
        touched.add(i)
touched_sat = touched - F0
print(f"\n97 handles touch {len(touched)} equations; {len(touched_sat)} currently-satisfied, {len(touched&F0)} failing(core)")
# distribution: how many handles per touched-sat eq
print("core eqs touched by handles:", sorted(touched & CORE))
