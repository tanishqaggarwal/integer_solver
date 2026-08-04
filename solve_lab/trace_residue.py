#!/usr/bin/env python3
import json
from collections import defaultdict
from propagate import NVARS
p=2**256-2**32-977
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
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
# trace x_20492 = x_16738 + x_36065; expand each
def show(t,depth=0,maxd=4):
    if depth>maxd: return
    g=gdef.get(t)
    val=vA[t]
    vs=f"{val}" if abs(val)<10**11 else f"BIG({len(str(abs(val)))}d)"
    print("  "*depth+f"x_{t} = {g[0][:55] if g else 'FREE/INPUT'}   [val={vs}]")
    if g:
        for u in g[1]:
            show(u,depth+1,maxd)
print("=== x_20492 tree (want x_20492 = CONST1 since x_19892=0) ===")
show(20492,0,3)
print(f"\nCONST1 = {CONST1}")
print(f"x_20492 = {vA[20492]}")
print(f"CONST1 - x_20492 = {CONST1-vA[20492]}  (%p = {(CONST1-vA[20492])%p})")
# check x_36065, x_16738 values vs known residues/consts
C1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
C2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
for nm,cv in [('C1',C1),('C2',C2),('CONST1',CONST1),('CONST2',CONST2)]:
    for t in [16738,36065,20492,10878,22542]:
        if vA[t]==cv or vA[t]%p==cv%p:
            print(f"  x_{t} matches {nm} (mod p: {vA[t]%p==cv%p}, exact: {vA[t]==cv})")
