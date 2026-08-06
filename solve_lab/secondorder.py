#!/usr/bin/env python3
"""Find product slacks x_a*x_b where BOTH factors are 0 & free (second-order DOF).
Focus on those in the obstruction cone (residue checks, verifiers, 16-ripple)."""
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
gate_out=set(); consumers=defaultdict(int)
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gate_out.add(d['t'])
        for u in d['vids']: consumers[u]+=1
freeinp=set(range(NVARS))-gate_out
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
vA={}
d=json.load(open('best_agentA_39022.json'))
for k,val in d.items(): vA[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
def val(x): return vA.get(x,0)
# atom membership
mem=defaultdict(int)
for pp in atoms:
    for x in atom_vars(pp): mem[x]+=1
# both-zero product slacks: monomial (a,b), val(a)==0 and val(b)==0
count=0; bothfree=0; clean=[]
for ai,poly in enumerate(atoms):
    for m in poly:
        if len(m)==2:
            a,b=m
            if a!=b and val(a)==0 and val(b)==0 and a not in wire and b not in wire:
                count+=1
                af=a in freeinp; bf=b in freeinp
                if af and bf:
                    bothfree+=1
                    # clean = both low atom-count (private-ish)
                    if mem[a]<=6 and mem[b]<=6:
                        clean.append((ai,a,b,mem[a],mem[b]))
print(f"both-zero non-wire product monomials: {count}")
print(f"  with BOTH factors free: {bothfree}")
print(f"  with both free & low-membership (<=6 atoms each): {len(clean)}")
# which touch the obstruction atoms/vars?
OBST={3277,3279,44601,20862,20864,42669,45276,44342,45677}
obstvars={6418,12553,2099,19964,4432,7068,642,28730,26777,13458}
touch=[c for c in clean if c[0] in OBST or (atom_vars(atoms[c[0]])&obstvars)]
print(f"\n  clean both-zero slacks in/near obstruction: {len(touch)}")
for ai,a,b,ma,mb in touch[:20]:
    print(f"    atom {ai}: x_{a}(mem{ma},free)*x_{b}(mem{mb},free) ; atom touches {sorted(atom_vars(atoms[ai])&obstvars)}")
