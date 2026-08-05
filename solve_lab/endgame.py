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
freeinp=set(range(NVARS))-set(gdef)
d=json.load(open('fc_partial.json'))
v=[0]*NVARS
for k,val in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
def tr(t,dep=0,md=3):
    g=gdef.get(t)
    print("  "*dep+f"x_{t}={g[0][:45] if g else 'FREE'} val={v[t] if abs(v[t])<10**8 else 'BIG'} free={t in freeinp} cons={len(consumers[t])}")
    if g and dep<md:
        for u in g[1]: tr(u,dep+1,md)
# divisibility check for 17897
print(f"atom 17897: need 13523997 | x_9106*x_21279. x_9106={v[9106]}, x_21279={v[21279]}")
print(f"  x_9106 % 13523997 = {v[9106]%13523997}, product%13523997 = {(v[9106]*v[21279])%13523997}")
print("\n--- x_9106 gate tree ---"); tr(9106,0,3)
print("\n--- x_31731 gate tree ---"); tr(31731,0,3)
print("\n--- x_9629 (sink) gate ---"); tr(9629,0,1)
