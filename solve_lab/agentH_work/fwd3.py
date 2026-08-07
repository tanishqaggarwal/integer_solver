"""Maximal orientation: use ANY atom with a unit-coefficient isolated variable as a definition."""
import model, pickle, os
from collections import Counter, defaultdict, deque
HERE=os.path.dirname(os.path.abspath(__file__))
d=model.get(); atom_src=d['atom_src']; atom_vars=d['atom_vars']
cands=pickle.load(open(os.path.join(HERE,'defcands.pkl'),'rb'))
NA=len(atom_src); NV=38748
iscand=defaultdict(list)   # var -> atoms where it is a def candidate
for a,cs in enumerate(cands):
    for v in cs: iscand[v].append(a)
prim=[v for v in range(NV) if not iscand[v]]
print('never a def candidate (primary free):', len(prim))
var_atoms=defaultdict(list)
for a,vs in enumerate(atom_vars):
    for v in vs: var_atoms[v].append(a)
known=bytearray(NV)
for v in prim: known[v]=1
definer=[-1]*NV; order=[]
# unknown-count per atom
unk=[0]*NA
for a in range(NA): unk[a]=sum(1 for v in atom_vars[a] if not known[v])
Q=deque(a for a in range(NA) if unk[a]==1)
while Q:
    a=Q.popleft()
    if unk[a]!=1: continue
    v=None
    for u in atom_vars[a]:
        if not known[u]: v=u; break
    if v is None or v not in cands[a]: continue
    known[v]=1; definer[v]=a; order.append(v)
    for b in var_atoms[v]:
        unk[b]-=1
        if unk[b]==1: Q.append(b)
nk=sum(known)
print('known:', nk, 'free inputs:', NV-nk)
used=set(definer[v] for v in range(NV) if definer[v]>=0)
print('checks:', NA-len(used))
