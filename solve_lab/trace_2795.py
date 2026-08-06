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
    dd=json.load(open(path)); v=[0]*NVARS
    for k,val in dd.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
freeinp=set(range(NVARS))-set(gdef)
def tr(t,d=0,maxd=4):
    g=gdef.get(t)
    val=vA[t]
    print("  "*d+f"x_{t}={g[0][:42] if g else 'FREE'} val={val if abs(val)<10**7 else 'BIG'} free={t in freeinp} consumers={len(consumers[t])}")
    if g and d<maxd:
        for u in g[1]: tr(u,d+1,maxd)
print("=== gap 24105 slack: x_25295 (want = 4261533*(x_31339-x_6858)) ===")
tr(25295)
print("\n=== gap 27902 slack: x_29967 (want = 12846437*(x_14853-x_1308)) ===")
tr(29967)
