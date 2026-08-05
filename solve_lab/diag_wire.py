#!/usr/bin/env python3
"""CRITICAL: is the wire (x_5101,x_32017,x_26789 in the core product-slacks) genuinely FORCED to p,
or free? Find single-variable atoms; build signed union-find over 2-term identity atoms; check what
forces x_26064/wire to p, and whether x_5101,x_32017,x_26789 are aliased into it."""
import json, re, ast, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
A=load_atoms()
# single-variable atoms: which variables are forced to a constant, and to what?
forced={}
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0:
            forced[v]=(-c0)//c1
print(f"single-var forced variables: {len(forced)}")
# how many forced to p?
fp=[v for v,x in forced.items() if x==p]
print(f"forced to EXACTLY p: {len(fp)} vars: {sorted(fp)[:10]}")
print(f"is x_26064 forced? {forced.get(26064)}  (==p: {forced.get(26064)==p})")
# are x_5101, x_32017, x_26789, x_4376, x_16123, x_35148 forced directly?
for v in [5101,32017,26789,4376,16123,35148,26064]:
    print(f"  x_{v}: forced={forced.get(v)} ({'=p' if forced.get(v)==p else forced.get(v)})")
# check the gate defs of x_5101 etc (are they identity-aliases?)
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gd={t:(rhs,vids) for t,rhs,vids in gates}
for v in [5101,32017,26789]:
    print(f"  gate x_{v} = {gd.get(v)}")
