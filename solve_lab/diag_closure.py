#!/usr/bin/env python3
"""Measure the handle-sharing closure of the 23 broken equations (how big is the coupled block)."""
import json, re, sys
from collections import defaultdict
from propagate import NVARS
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
eqfree=[vs&freeinp for vs in eqvars]
# map free input -> equations
byh=defaultdict(set)
for i,vs in enumerate(eqfree):
    for v in vs: byh[v].add(i)
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
seed=[3408, 3841, 4134, 4526, 5069, 7276, 15440, 15724, 15927, 21600, 22139, 22825, 27289, 27999, 28718, 29305, 31134, 31269, 32463, 33195, 36387, 36390, 38888]
# 1-hop and full closure over shared handles (excluding controls/quotients which we treat specially)
CTRL={14853,12186,16742,19750,30317,2936,5146}
eqs=set(seed); handles=set()
for _ in range(20):
    newh=set()
    for i in eqs: newh|=(eqfree[i]-CTRL)
    neweq=set()
    for h in newh: neweq|=byh[h]
    neweq-=CORE
    if neweq<=eqs and newh<=handles: break
    eqs|=neweq; handles|=newh
print(f"closure over shared handles: {len(eqs)} equations, {len(handles)} handles")
# 1-hop only
h1=set()
for i in seed: h1|=(eqfree[i]-CTRL)
e1=set()
for h in h1: e1|=byh[h]
e1-=CORE
print(f"1-hop: {len(e1)} equations, {len(h1)} handles")
