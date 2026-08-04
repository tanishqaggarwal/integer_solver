#!/usr/bin/env python3
"""Get x_24908 free-input cone size, and check overlap of critical-pin set with the 23 broken eqs' handles."""
import json, re, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
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
fc=freecone(24908)
print(f"x_24908 free cone: {len(fc)} free inputs")
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
broken=[3408, 3841, 4134, 4526, 5069, 7276, 15440, 15724, 15927, 21600, 22139, 22825, 27289, 27999, 28718, 29305, 31134, 31269, 32463, 33195, 36387, 36390, 38888]
crit={14853,12186,16742}|fc
bh=set()
for i in broken: bh|=(eqvars[i]&freeinp)
print(f"broken-eq free handles: {len(bh)}; overlap with critical pin set: {len(bh&crit)}")
print(f"  overlap = {sorted(bh&crit)}")
# also does x_24908 cone overlap the broken handles a lot?
print(f"x_24908 cone ∩ broken handles: {sorted(fc&bh)[:20]}")
