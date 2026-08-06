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
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.val[13195]=1; H.forward()
def status(v):
    isfree=v in H.freeinp; rhs=''
    if v in H.definer: rhs=H.gates[H.definer[v]][1][:45]
    return ('FREE' if isfree else 'gate'),len(var_in_atom.get(v,[])),rhs,H.val[v]%p
# pin targets and gadget vars
print("=== pin targets & gadget vars (x_13195=1) ===")
for v in [18623,6467,30175,37405,14608,33441,19801,1279,35654,19247,6418,12553,2803,3082]:
    st,na,rhs,val=status(v)
    print(f"x_{v}: {st} #at={na} val={'0' if val==0 else 'nz'} rhs={rhs}")
# free ancestors of gadget product vars
for v in [30175,19247,37405,14608,33441]:
    anc=H.anc.get(v,{v})
    fa=[a for a in anc if a in H.freeinp]
    print(f"x_{v}: {len(fa)} free ancestors; can set to 0? (self free={v in H.freeinp})")
