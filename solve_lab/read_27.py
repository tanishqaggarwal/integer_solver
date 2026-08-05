#!/usr/bin/env python3
import heal_harness as H
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=H.p
atoms=load_atoms()
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
mem=defaultdict(int)
for pp in atoms:
    for x in atom_vars(pp): mem[x]+=1
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=v013.get(v,0)
H.val[14853]=vA[14853]; H.val[31339]=vA[31339]
H.forward(); v=H.val
# gap atoms 24105 (x_14853... wait it said 31339), 27902
for ai in [24105,27902]:
    print(f"=== atom {ai} (gap) ===")
    for m,c in sorted(atoms[ai].items(),key=lambda kv:(len(kv[0]),kv[0])):
        vs="*".join(f"x_{x}" for x in m) if m else "1"
        meta=[f"x{x}[free={x in freeinp},atoms={mem[x]},wire={x in wire},v={v[x] if abs(v[x])<10**9 else 'BIG'}]" for x in m]
        print(f"   {c:+d} * {vs}  {meta}")
# check the fine-grained slack factors
print("\n=== fine-grained slack factor values ===")
for zf,nf in [(8183,12378),(30060,12378),(14515,36977),(19750,36977),(5616,30033),(19275,15658)]:
    print(f"  x_{zf}(free={zf in freeinp},v={v[zf]}) * x_{nf}(free={nf in freeinp},wire={nf in wire},v={v[nf] if abs(v[nf])<10**9 else 'BIG'}, v%p={v[nf]%p if abs(v[nf])>10**9 else '-'})")
