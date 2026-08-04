#!/usr/bin/env python3
import json,re
from propagate import NVARS
p=2**256-2**32-977
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
VAR=re.compile(r'x_(\d+)')
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
def fails(v):
    ns={'v':v,'__builtins__':{}}
    return [i for i,c in enumerate(eqcode) if eval(c,ns)!=0]
base=set(fails(vA))
print(f"agentA base: {len(base)} fails")
# Test G2 split: x_17499 -> 1, x_9413 = x_4432 - x_19964, x_28730 = that
def test(name, changes):
    v=vA[:]
    for idx,val in changes: v[idx]=val
    F=set(fails(v))
    print(f"{name}: {len(F)} fails (fixed {len(base-F)}: {sorted(base-F)}, broke {len(F-base)}: {sorted(F-base)[:20]})")
    return F
g2=[(17499,1),(9413,vA[4432]-vA[19964]),(28730,vA[4432]-vA[19964])]
test("G2 split (x_17499=1)", g2)
# G1 split: x_28599 -> 1, x_642 = x_17325, need 7376877 | (x_7068-x_2099)
diff=vA[7068]-vA[2099]
print(f"\nG1: (x_7068-x_2099) mod 7376877 = {diff%7376877} (need 0 for integer x_17325)")
# try both splits together assuming G1 divisibility handled by x_17325 (may be non-integer -> skip if not)
if diff%7376877==0:
    x17325=diff//7376877
    g1=[(28599,1),(17325,x17325),(642,x17325)]
    test("G1 split (x_28599=1)", g1)
    test("BOTH splits", g2+g1)
else:
    print("G1 needs x_7068 or x_2099 adjusted for divisibility; testing G2 split alone above.")
