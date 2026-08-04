#!/usr/bin/env python3
"""Find defs of x_3558, x_29322. Trace to settable free partner whose residue controls them linearly."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
def show(v,depth=0,seen=None):
    if seen is None: seen=set()
    ind='  '*depth
    if v in freeinp:
        print(f"{ind}x_{v} [FREE]"); return
    if v in seen or depth>4:
        print(f"{ind}x_{v} = {gate_defs[v][0][:50]} ..."); return
    seen.add(v)
    rhs,vids=gate_defs[v]
    print(f"{ind}x_{v} = {rhs[:60]}")
    for u in vids: show(u,depth+1,seen)
for root in [3558,29322]:
    print(f"\n===== x_{root} definition tree =====")
    show(root)
