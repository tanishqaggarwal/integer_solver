import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
# atom fanout: which atoms contain each var
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
var_in_atom=defaultdict(list)
for ai,a in enumerate(atoms):
    vs=set()
    for mono,c in a['poly']:
        for v in mono: vs.add(v)
    for v in vs: var_in_atom[v].append(ai)

kv=[7068,2099,642,17325,4432,19964,28730,9413,6418,12553,31861,14865,9118,8731,2081,4287]
print("var | free? | #atoms | gate_rhs(if gate)")
gdef={t:(rhs,vids) for t,rhs,vids in H.gates}
# definer maps target->gate index in H.gates
for v in kv:
    isfree = v in H.freeinp
    natoms=len(var_in_atom.get(v,[]))
    rhs=''
    if v in H.definer:
        gi=H.definer[v]; t,r,vids=H.gates[gi]; rhs=r[:60]
    print(f"x_{v} | {'FREE' if isfree else 'gate'} | {natoms} | {rhs}")

print("\n--- atoms containing 7068 ---")
for ai in var_in_atom[7068]:
    print(ai, atoms[ai]['repr'][:90], " eqs=",atoms[ai]['eqs'][:6])
print("\n--- atoms containing 4432 ---")
for ai in var_in_atom[4432]:
    print(ai, atoms[ai]['repr'][:90], " eqs=",atoms[ai]['eqs'][:6])
