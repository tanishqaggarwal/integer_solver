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
# full repr of core verifier atoms
for ai in [33109,32180,44619,15517,15518,15521,15523,15516,15520]:
    print(f"--- atom {ai} (eqs {sorted(set(atoms[ai]['eqs']))[:4]}) ---")
    print("  ", atoms[ai]['repr'][:300])
# status of core-feeding vars
def status(v):
    isfree=v in H.freeinp
    rhs=''
    if v in H.definer:
        gi=H.definer[v]; t,r,vids=H.gates[gi]; rhs=r[:45]
    return ('FREE' if isfree else 'gate'), rhs
print("\n=== core-feeding vars ===")
for v in [9192,24453,30213,22162,17601,29356,17702,27762,11602,32680,5647,15298]:
    st,rhs=status(v); print(f"x_{v}: {st} rhs={rhs}")
# QR check at agentA
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
for v in [33469,29322,3558,27713,1326]:
    print(f"x_{v} mod p = {H.val[v]%p}")
x=H.val[33469]%p
qr=pow(x,(p-1)//2,p)
print(f"x_33469 QR? legendre={qr} ({'QR' if qr==1 else 'NQR' if qr==p-1 else 'zero'})")
