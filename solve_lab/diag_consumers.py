#!/usr/bin/env python3
"""What consumes x_24908, x_29322, x_16742, x_14853, x_12186? Find their forward-cone of gate
consumers and whether those feed NON-core equations (coupling) or only core (private)."""
import json, re, sys
from collections import defaultdict, deque
from propagate import NVARS
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}; consumers=defaultdict(list)
for t,rhs,vids in gates:
    gate_defs[t]=(rhs,vids)
    for u in vids: consumers[u].append(t)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
byvar=defaultdict(set)
for i,vs in enumerate(eqvars):
    for v in vs: byvar[v].add(i)
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
def fwd_gate_cone(root):
    seen=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        for t in consumers.get(x,()): st.append(t)
    return seen
for root in [24908, 29322, 16742, 14853, 12186, 3558]:
    fc=fwd_gate_cone(root)  # all gates (and root) reachable forward
    eqs=set()
    for v in fc: eqs|=byvar.get(v,set())
    noncore=sorted(eqs-CORE); core=sorted(eqs&CORE)
    print(f"x_{root}: forward-reaches {len(fc)} vars; appears(via cone) in {len(core)} core + {len(noncore)} non-core eqs")
    if len(noncore)<=25: print(f"    non-core: {noncore}")
