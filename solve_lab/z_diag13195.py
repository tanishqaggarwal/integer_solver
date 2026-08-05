import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
eq_atoms=defaultdict(list)
for ai,a in enumerate(atoms):
    for e in set(a['eqs']): eq_atoms[e].append(ai)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.val[13195]=1
H.forward()
ns={'v':H.val,'__builtins__':{}}
Fp=[i for i,c in enumerate(H.eqcode) if eval(c,ns)%p!=0]
print(f"x_13195=1: {len(Fp)} mod-p fails")
# which atoms are nonzero across these fails
def ev_atom(a):
    s=0
    for mono,c in a['poly']:
        t=c
        for v in mono: t*=H.val[v]
        s+=t
    return s%p
nz=set()
for e in Fp:
    for ai in eq_atoms[e]:
        if ev_atom(atoms[ai])!=0: nz.add(ai)
print(f"nonzero atoms: {len(nz)}")
for ai in sorted(nz):
    print(f"  atom {ai}: {atoms[ai]['repr'][:70]}")
# loads
for nm,v in [('L1',11150),('L2',25739),('L3',37758)]:
    print(f"{nm}=x_{v} mod p = {H.val[v]%p}")
