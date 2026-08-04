#!/usr/bin/env python3
"""Fix G2 via the FINE-GRAINED slack x_19892 = x_21279*x_8731 (x_21279=x_9062*x_20434),
making x_19964 = x_4432 exactly. Test ripple."""
import json
from collections import defaultdict
from propagate import NVARS
p=2**256-2**32-977
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
for lbl,v in [('x_9062',9062),('x_20434',20434),('x_8731',8731),('x_21279',21279),('x_19892',19892),('x_20492',20492),('x_19964',19964),('x_4432',4432)]:
    print(f"  {lbl}: free={v in freeinp}, val={vA[v] if abs(vA[v])<10**11 else 'BIG(%d dig)'%len(str(abs(vA[v])))}")
# The gap to absorb:
G = vA[4432]-vA[20492]   # since x_19964=x_20492 currently, want x_19964=x_4432 => x_19892 = x_4432-x_20492
print(f"\ngap x_4432 - x_20492 = {G}")
print(f"  need x_19892 = G, via x_21279*x_8731. Set x_8731=1, x_21279=G, via x_9062=G (x_20434=1).")
# check x_9062, x_8731, x_20434 free
print(f"  x_9062 free={9062 in freeinp}, x_8731 free={8731 in freeinp}, x_20434 free={20434 in freeinp}")
# Load equations
import re
VAR=re.compile(r'x_(\d+)')
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
def fails(v):
    ns={'v':v,'__builtins__':{}}
    return set(i for i,c in enumerate(eqcode) if eval(c,ns)!=0)
base=fails(vA)
print(f"\nagentA base fails: {len(base)}")
# Construct: set free inputs x_9062=G, x_8731=1 (x_20434 already 1), then set the gate chain consistently
v=vA[:]
v[9062]=G       # if free
v[8731]=1
# recompute the chain manually (and everything downstream would need forward; but test direct chain first)
v[21279]=v[9062]*v[20434]   # = G*1 = G
v[19892]=v[21279]*v[8731]   # = G*1 = G
v[19964]=v[20492]+v[19892]  # = x_20492 + G = x_4432
F=fails(v)
print(f"after chain set (direct, no full forward): {len(F)} fails")
print(f"  fixed: {sorted(base-F)}")
print(f"  broke: {sorted(F-base)[:30]}")
