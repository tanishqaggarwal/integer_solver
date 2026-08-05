#!/usr/bin/env python3
import heal_harness as H
import json
from collections import defaultdict
p=H.p
gdef={}; consumers=defaultdict(list)
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gdef[d['t']]=(d['rhs'],tuple(d['vids']))
        for u in d['vids']: consumers[u].append(d['t'])
vA=H.loadd('best_agentA_39022.json')
R12553_old=42775533402728869434716629464193396056515231264222641773817154079369026410240838606908039
R6418_old=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
gap12553=CONST1-R12553_old
g1=CONST2-R6418_old
print(f"atom 3279 wants x_13458 = {gap12553}")
print(f"atom 3277 wants 15804267*x_26777 = {g1}; g1 % 15804267 = {g1%15804267} (divisible: {g1%15804267==0})")
# trace x_13458, x_26777 gates
def tr(t,d=0):
    g=gdef.get(t)
    print("  "*d+f"x_{t} = {g[0][:50] if g else 'FREE'}  free={t not in gdef} consumers={len(consumers[t])}")
    if g and d<3:
        for u in g[1]: tr(u,d+1)
print("\n=== x_13458 gate ==="); tr(13458)
print("\n=== x_26777 gate ==="); tr(26777)
