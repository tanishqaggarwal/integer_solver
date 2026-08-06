import heal_harness as H
from propagate import load_atoms, atom_vars
from collections import defaultdict, deque
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V=H.val
A=load_atoms()
var_atoms=defaultdict(list)
for ai,poly in enumerate(A):
    for v in atom_vars(poly): var_atoms[v].append(ai)
def deg(poly): return max(len(m) for m in poly)
def evalp(poly):
    s=0
    for m,c in poly.items():
        t=c
        for v in m: t*=V[v]
        s+=t
    return s
# The wire class (vars ~ p) and their identification
wire=set(v for v in range(len(V)) if V[v]!=0 and V[v]%p==0)
print("wire-ish vars (val = mult of p, nonzero):", len(wire))
# BFS closure: start from vars in G1,G2 that are NOT slack multipliers / not wire.
# A var is 'movable'. Moving it breaks all atoms containing it. To heal each such atom,
# we need to move ANOTHER var in it. Track reachable atoms & vars.
# We treat wire vars and rare-partner slacks specially (terminators).
seed_vars={4432,7068,2099,19964,642,28730}  # obstruction vars
# free inputs are the true DOF; gates are determined. But atoms constrain everything.
# Do closure over NON-wire, NON-const vars.
def is_frozen(v):
    # vars we won't move: wire members (val=mult of p), and huge structural consts
    return v in wire
seen_atoms=set()
seen_vars=set()
q=deque(seed_vars)
DEGS=defaultdict(int)
while q:
    v=q.popleft()
    if v in seen_vars: continue
    seen_vars.add(v)
    for ai in var_atoms[v]:
        if ai in seen_atoms: continue
        seen_atoms.add(ai)
        poly=A[ai]; d=deg(poly); DEGS[d]+=1
        for u in atom_vars(poly):
            if u not in seen_vars and not is_frozen(u):
                q.append(u)
print(f"Closure: {len(seen_atoms)} atoms, {len(seen_vars)} vars")
print("atom degrees in closure:", dict(DEGS))
# How many of these atoms are currently NONZERO?
nz=[ai for ai in seen_atoms if evalp(A[ai])!=0]
print("nonzero atoms in closure:", nz)
# free inputs in closure
fi=[v for v in seen_vars if v in H.freeinp]
print(f"free inputs in closure: {len(fi)}")
