import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
var_in_atom=defaultdict(list)
for ai,a in enumerate(atoms):
    vs=set(v for mono,c in a['poly'] for v in mono)
    for v in vs: var_in_atom[v].append(ai)
def status(v):
    isfree=v in H.freeinp
    rhs=''
    if v in H.definer:
        gi=H.definer[v]; t,r,vids=H.gates[gi]; rhs=r[:55]
    return ('FREE' if isfree else 'gate'), len(var_in_atom.get(v,[])), rhs
print("=== CORE vars ===")
core=[33469,29322,3558,27713,1326,14853,12186,24908,16742]
for v in core:
    st,na,rhs=status(v); print(f"x_{v}: {st} #atoms={na} rhs={rhs}")
print("\n=== atoms mentioning core S,T (33469,29322,3558,27713) ===")
for v in [33469,3558,27713]:
    print(f"-- x_{v} atoms --")
    for ai in var_in_atom[v]:
        print(f"  {ai}: {atoms[ai]['repr'][:75]}")
