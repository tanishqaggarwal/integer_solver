#!/usr/bin/env python3
import heal_harness as H
import json
from collections import defaultdict
p=H.p
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gdef[d['t']]=(d['rhs'],tuple(d['vids']))
vA=H.loadd('best_agentA_39022.json')
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
# trace x_2099 = x_37158 + x_25297 (x_25297=0), x_37158 = x_10878 + x_22542
def leafres(t,d=0):
    g=gdef.get(t)
    if not g or d>4: return t
    print("  "*d+f"x_{t}={g[0][:45]} val={str(vA.get(t,0))[:8] if abs(vA.get(t,0))<10**11 else 'BIG'} free={t in H.freeinp}")
    for u in g[1]: leafres(u,d+1)
print("=== x_2099 tree ===")
leafres(2099)
# find which free input equals x_2099's value (the residue holding x_2099)
print(f"\nx_2099 val = {str(vA[2099])[:12]}..")
# candidates: free inputs whose value == x_2099
cand=[v for v in H.freeinp if vA.get(v,0)==vA[2099]]
print(f"free inputs holding x_2099's residue: {cand}")
cand2=[v for v in H.freeinp if vA.get(v,0)==vA[19964]]
print(f"free inputs holding x_19964's residue: {cand2}")
