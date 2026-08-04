#!/usr/bin/env python3
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
# atom membership count per var
mem=defaultdict(list)
for ai,poly in enumerate(atoms):
    for v in atom_vars(poly): mem[v].append(ai)
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json'); v013=loadv('best/new_instance_partial_39013.json')
# gate outputs / free inputs
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'],d['rhs'],tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
interest={'x_642(slack1)':642,'x_28599(fac)':28599,'x_17325(fac)':17325,
          'x_28730(slack2)':28730,'x_17499(fac)':17499,'x_9413(fac)':9413,
          'x_2099':2099,'x_7068':7068,'x_4432':4432,'x_19964':19964}
for lbl,v in interest.items():
    print(f"{lbl}: free={v in freeinp}, in {len(mem[v])} atoms, valA={vA[v]}, val013={v013[v]}")
print("\n=== divisibility for G1: need 7376877 | (x_7068 - x_2099) ===")
diff=vA[7068]-vA[2099]
print(f"x_7068 - x_2099 (agentA) = {diff}")
print(f"  mod 7376877 = {diff%7376877}  (divisible: {diff%7376877==0})")
diff013=v013[7068]-v013[2099]
print(f"x_7068 - x_2099 (39013) = {diff013}, mod 7376877 = {diff013%7376877}")
print(f"x_642 at 39013 = {v013[642]}, 7376877*x_642^013 = {7376877*v013[642]}")
print(f"\n=== G2 (coeff 1, trivial): x_28730 target = x_4432 - x_19964 (agentA) = {vA[4432]-vA[19964]}")
# factor 7376877
n=7376877
f=[]; d=2; m=n
while d*d<=m:
    while m%d==0: f.append(d); m//=d
    d+=1
if m>1: f.append(m)
print(f"7376877 factors = {f}")
