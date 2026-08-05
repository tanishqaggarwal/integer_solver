#!/usr/bin/env python3
"""Examine the 13 equations broken by the x_16742 fix. Are they definitions, load-packings, or
verifier squares? What loads do they contain and what's their residual structure?"""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
LOADS={11150,25739,37758}
broken=[3408, 3841, 4134, 4526, 5069]
for i in broken:
    L=lines[i].rsplit('=',1)[0]
    vs=set(int(m) for m in VAR.findall(L))
    loads_in=vs&LOADS
    # is it a square? outermost is c*(root); check root is X*X
    print(f"\n=== eq {i}: len={len(L)}, loads={loads_in}, #vars={len(vs)}")
    # detect leading structure
    print(f"    head: {L[:230]}")
