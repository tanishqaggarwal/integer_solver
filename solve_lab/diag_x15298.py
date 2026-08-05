#!/usr/bin/env python3
"""How is x_15298 determined? What sets it to 1? Check x_7715, x_34554 and whether forcing
x_15298=0 is possible (any equation forcing it to 1?). Count core eqs that vanish if x_15298=0."""
import json, re, ast, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
print("x_15298 def:", gate_defs.get(15298))
for v in [7715,34554]:
    print(f"x_{v}: {'FREE' if v in freeinp else gate_defs.get(v)}")
# what feeds x_7715, x_34554 (backward)
def show(v,d=0,seen=None):
    if seen is None: seen=set()
    if v in freeinp: print("  "*d+f"x_{v} FREE"); return
    if v in seen or d>3: print("  "*d+f"x_{v}=..."); return
    seen.add(v); print("  "*d+f"x_{v} = {gate_defs[v][0][:50]}")
    for u in gate_defs[v][1]: show(u,d+1,seen)
print("--- x_7715 tree ---"); show(7715)
print("--- x_34554 tree ---"); show(34554)
# is x_15298 (or x_7715/x_34554) directly constrained by an equation like x^2=x (boolean) or =1?
occ15298=[i for i in range(len(lines)) if 15298 in set(int(m) for m in VAR.findall(lines[i]))]
print(f"x_15298 appears in {len(occ15298)} equations")
# count how many core-eq terms are x_15298-multiplied
CORE=[2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]
print(f"x_15298 in core eqs: {[i for i in CORE if i in occ15298]}")
