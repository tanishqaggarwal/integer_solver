import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
# map eq -> atoms
eq_atoms=defaultdict(list)
for ai,a in enumerate(atoms):
    for e in set(a['eqs']): eq_atoms[e].append(ai)

vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
x2099=H.val[2099]; x19964=H.val[19964]
H.val[7068]=x2099; H.val[4432]=x19964
H.forward()

def ev_atom(a):
    s=0
    for mono,c in a['poly']:
        t=c
        for v in mono: t*=H.val[v]
        s+=t
    return s

broken=[697, 1985, 5225, 10815, 16048, 17784, 17801, 22402, 23667, 24721, 27124, 28737, 29638, 29959, 35935, 37431]
nz_atoms=set()
for e in broken:
    for ai in eq_atoms[e]:
        val=ev_atom(atoms[ai])
        if val%p!=0 or val!=0:
            nz_atoms.add(ai)
print(f"nonzero atoms across 16 broken eqs: {len(nz_atoms)}")
for ai in sorted(nz_atoms):
    val=ev_atom(atoms[ai])
    print(f"atom {ai}: {atoms[ai]['repr'][:80]}")
    print(f"      value mod p = {val%p}   ==0inZ:{val==0}  eqs={sorted(set(atoms[ai]['eqs']))[:8]}")
