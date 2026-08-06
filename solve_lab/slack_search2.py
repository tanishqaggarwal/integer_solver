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
vA=H.loadd('best_agentA_39022.json')
def val(x): return vA.get(x,0)
# For x_6418 and x_12553: every product slack (a*b) in any atom containing them, with factor values
for tgt in [6418,12553]:
    print(f"\n=== residue x_{tgt}: product-slack factors across all its atoms ===")
    seen=set()
    for ai,poly in enumerate(atoms):
        if tgt not in atom_vars(poly): continue
        for m in poly:
            if len(m)==2 and tgt not in m:  # a product NOT involving tgt (a potential compensating slack in same atom)
                a,b=m
                key=tuple(sorted(m))
                if key in seen: continue
                seen.add(key)
                # is it activatable & fine-grained? one factor free&0, other's value
                for zf,nf in [(a,b),(b,a)]:
                    if zf in freeinp and val(zf)==0:
                        nfv=val(nf); nfwire = nf in wire
                        gran = 'p-gran' if (nfv%p==0 and nfv!=0) or nfwire else ('FINE' if nfv!=0 else 'zero')
                        if gran=='FINE':
                            print(f"  atom {ai}: FINE slack x_{zf}(free,0) * x_{nf}(val={nfv if abs(nfv)<10**9 else 'BIG'}) [gran={gran}]")
    # also: slacks that DIRECTLY set tgt-related residue-load (coeff on a private handle)
print("\n=== is CONST1/CONST2 a loadable HUGE constant? (search equation text) ===")
import subprocess
for nm,c in [('CONST1',97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680),
             ('CONST2',126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506)]:
    n=sum(1 for L in open('../EQUATIONS.txt') if str(c) in L)
    print(f"  {nm} appears literally in {n} equations")
