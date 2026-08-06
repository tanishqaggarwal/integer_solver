#!/usr/bin/env python3
import json
from collections import defaultdict
from propagate import NVARS
p=2**256-2**32-977
gdef={}; consumers=defaultdict(list)
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gdef[d['t']]=(d['rhs'],tuple(d['vids']))
        for u in d['vids']: consumers[u].append(d['t'])
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
# Does x_2081 feed x_15298 (the MUX)? trace x_15298 ancestors
def ancestors(t, seen=None):
    if seen is None: seen=set()
    if t in seen: return seen
    seen.add(t)
    if t in gdef:
        for u in gdef[t][1]: ancestors(u,seen)
    return seen
anc15298=ancestors(15298)
print(f"x_2081 in x_15298 ancestors (MUX)? {2081 in anc15298}")
print(f"x_24601 in x_15298 ancestors? {24601 in anc15298}")
print(f"x_15298 free ancestors involving activators: {[a for a in anc15298 if a in (2081,24601,8599,21839,25956,7304)]}")
print(f"x_2081 direct consumers: {consumers[2081]}")
for c in consumers[2081][:12]:
    print(f"   x_{c} = {gdef[c][0][:45]}")
print(f"\nx_2081 val={vA[2081]}, x_24601 val={vA[24601]}, x_15298 val={vA[15298]}")
# What feeds x_15298?
print(f"x_15298 = {gdef.get(15298,('?',()))[0]}")
for u in gdef.get(15298,('',()))[1]:
    print(f"   x_{u} = {gdef.get(u,('FREE',()))[0][:40]}")
