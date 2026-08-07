"""Canonical forward orientation: greedily orient atoms into definitions; find free inputs."""
import model, pickle, os, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
d = model.get()
atom_src=d['atom_src']; atom_vars=d['atom_vars']; eq_terms=d['eq_terms']
polys = pickle.load(open(os.path.join(HERE,'polys.pkl'),'rb'))
cands = pickle.load(open(os.path.join(HERE,'defcands.pkl'),'rb'))
NA=len(atom_src); NV=38748

# For each atom, list of vars; for each var, atoms containing it
var_atoms = defaultdict(list)
for a,vs in enumerate(atom_vars):
    for v in vs: var_atoms[v].append(a)

# greedy triangularization: an atom is "ready" if exactly one of its vars is undetermined
# and that var is a def candidate.
undet = [len(vs) for vs in atom_vars]   # number of undetermined vars in atom
determined = bytearray(NV)
definer = [-1]*NV      # atom that defines var
order = []             # topological order of definitions
used_as_def = bytearray(NA)

from collections import deque
# ready queue: atoms with undet count == 1 whose remaining var is a def candidate
def remaining_var(a):
    for v in atom_vars[a]:
        if not determined[v]: return v
    return None

Q = deque(a for a in range(NA) if undet[a]==1)
# also atoms with undet==0 are checks (fully determined) - nothing to do
while Q:
    a = Q.popleft()
    if undet[a]!=1: continue
    v = remaining_var(a)
    if v is None: continue
    if v not in cands[a]: continue
    if determined[v]: continue
    determined[v]=1; definer[v]=a; used_as_def[a]=1; order.append(v)
    for b in var_atoms[v]:
        undet[b]-=1
        if undet[b]==1: Q.append(b)

ndef = sum(determined)
print('defined vars:', ndef, 'free inputs:', NV-ndef)
free = [v for v in range(NV) if not determined[v]]
checks = [a for a in range(NA) if not used_as_def[a]]
print('check atoms (not used as definition):', len(checks))
print('atoms still with undetermined vars:', sum(1 for a in range(NA) if undet[a]>0))
pickle.dump({'determined':bytes(determined),'definer':definer,'order':order,
             'free':free,'checks':checks,'used_as_def':bytes(used_as_def)},
            open(os.path.join(HERE,'fwd.pkl'),'wb'))
