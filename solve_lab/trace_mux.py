#!/usr/bin/env python3
import json
from propagate import NVARS
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gdef[d['t']]=(d['rhs'],tuple(d['vids']))
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
def tr(t,d=0,maxd=5):
    g=gdef.get(t)
    print("  "*d+f"x_{t}={g[0][:40] if g else 'FREE'} val={vA[t] if abs(vA[t])<10**6 else 'BIG'}")
    if g and d<maxd:
        for u in g[1]: tr(u,d+1,maxd)
print("=== x_15298 = x_7715 * x_34554 -> need both=1 ===")
print("--- x_7715 branch ---"); tr(7715,0,5)
print("--- x_34554 branch ---"); tr(34554,0,5)
