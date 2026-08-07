"""Syntactic-target gate DAG.  target(atom) = leading bare variable of `X - rest`."""
import model, pickle, os, ast
from collections import Counter, defaultdict, deque
HERE=os.path.dirname(os.path.abspath(__file__))
d=model.get(); atom_src=d['atom_src']; atom_vars=d['atom_vars']; eq_terms=d['eq_terms']
NA=len(atom_src); NV=38748

def target_of(s):
    """atom source 'x_123 - rest'  -> 123 ; else None"""
    t = ast.parse(s, mode='eval').body
    # strip leading unary/leading const-mult? no: only pure top-level Sub with Name left
    if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Sub) and isinstance(t.left, ast.Name):
        return int(t.left.id[2:])
    return None

tgt = [target_of(s) for s in atom_src]
print('atoms with syntactic target:', sum(1 for t in tgt if t is not None))
istarget = set(t for t in tgt if t is not None)
print('distinct target vars:', len(istarget))
free0 = [v for v in range(NV) if v not in istarget]
print('never a target (pure free inputs):', len(free0))

var_atoms = defaultdict(list)
for a,vs in enumerate(atom_vars):
    for v in vs: var_atoms[v].append(a)

known = bytearray(NV)
for v in free0: known[v]=1
definer=[-1]*NV; order=[]
# count unknown non-target vars per atom
unk = [0]*NA
for a in range(NA):
    t=tgt[a]
    unk[a]=sum(1 for v in atom_vars[a] if v!=t and not known[v])
Q=deque(a for a in range(NA) if tgt[a] is not None and unk[a]==0)
while Q:
    a=Q.popleft()
    t=tgt[a]
    if t is None or known[t] or unk[a]!=0: continue
    known[t]=1; definer[t]=a; order.append(t)
    for b in var_atoms[t]:
        if tgt[b]!=t:
            unk[b]-=1
            if unk[b]==0 and tgt[b] is not None and not known[tgt[b]]: Q.append(b)
nk=sum(known)
print('known after propagation:', nk, ' unresolved (cycles):', NV-nk)
used=set(definer[v] for v in range(NV) if definer[v]>=0)
checks=[a for a in range(NA) if a not in used]
print('definitions used:', len(used), 'check atoms:', len(checks))
pickle.dump({'tgt':tgt,'free0':free0,'definer':definer,'order':order,
             'known':bytes(known),'checks':checks},open(os.path.join(HERE,'fwd2.pkl'),'wb'))
# check-atom degree/eq stats
eq_of_atom=defaultdict(list)
for i,(m,sq,tl) in enumerate(eq_terms):
    for c,a in tl: eq_of_atom[a].append(i)
print('checks in >1 eq:', sum(1 for a in checks if len(eq_of_atom[a])>1))
print('mean eqs per check:', sum(len(eq_of_atom[a]) for a in checks)/len(checks))
