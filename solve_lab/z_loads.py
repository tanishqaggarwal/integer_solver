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
    for v in set(v for mono,c in a['poly'] for v in mono): var_in_atom[v].append(ai)
pins=json.load(open('pinrec.json'))
selectors=set(r[1] for r in pins)
# loads
loads={'L1':11150,'L2':25739,'L3':37758}
for nm,v in loads.items():
    print(f"=== {nm} = x_{v}  free={v in H.freeinp} #atoms={len(var_in_atom[v])} ===")
    if v in H.definer:
        gi=H.definer[v]; t,r,vids=H.gates[gi]; print(f"  gate rhs: {r[:80]}")
    # bits in ancestry
    anc=H.anc.get(v,{v}); bits=anc&selectors
    print(f"  free ancestors: {len(anc)}, bits: {sorted(bits)}")
    for ai in var_in_atom[v][:6]:
        print(f"  atom {ai}: {atoms[ai]['repr'][:75]}")
# check load values at agentA
vA=H.loadd('best_agentA_39022.json')
for x in H.freeinp: H.val[x]=vA.get(x,0)
H.forward()
for nm,v in loads.items():
    print(f"{nm}=x_{v} at agentA mod p = {H.val[v]%p}, mod 6672769 = {H.val[v]%6672769}")
