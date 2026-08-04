#!/usr/bin/env python3
"""Measure the coupled component around G1/G2 in the forward model.
A check atom is 'active' if, at agentA, perturbing a free input in the cone changes it.
We approximate the cone by graph reachability using free-input ancestors (anc) of gate vars."""
import heal_harness as H
import sg2_lib as L
import pickle
from collections import defaultdict, deque
p = H.p
atoms = L.load_atoms_full()
A = {a['idx']: a for a in atoms}
idx = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/idx.pkl','rb'))
var_atoms = idx['var_atoms']

anc = H.anc  # gate var -> set of free-input ancestors
freeinp = H.freeinp

def free_deps(v):
    if v in freeinp: return {v}
    return anc.get(v, set())

# free-input deps of an atom
def atom_free_deps(ai):
    s = set()
    for v in L.atom_vars(A[ai]['poly']):
        s |= free_deps(v)
    return s

# For each free input, which atoms depend on it (contain a var whose free-anc includes it)
# Build: free input -> set of atoms.  (Precompute atom -> free deps once.)
atom_fdeps = {}
for a in atoms:
    s = set()
    for v in L.atom_vars(a['poly']):
        s |= free_deps(v)
    atom_fdeps[a['idx']] = s
fin_to_atoms = defaultdict(set)
for ai, s in atom_fdeps.items():
    for f in s:
        fin_to_atoms[f].add(ai)

# seed: free deps of G1 (atom 20862) and G2 (atom 20864)
seedatoms = [20862, 20864]
S = set()
for ai in seedatoms:
    S |= atom_fdeps[ai]
print(f"seed free inputs (from G1,G2): {len(S)}")

# closure: BFS over free inputs coupled via shared atoms
# Two free inputs are coupled if some atom depends on both.
active_atoms = set()
q = deque(S)
seen = set(S)
while q:
    f = q.popleft()
    for ai in fin_to_atoms[f]:
        if ai in active_atoms: continue
        active_atoms.add(ai)
        for g in atom_fdeps[ai]:
            if g not in seen:
                seen.add(g); q.append(g)

print(f"coupled free inputs: {len(seen)}")
print(f"active atoms (touching the cone): {len(active_atoms)}")
# how many are deg-1 vs deg-2 vs deg-4
from collections import Counter
degc = Counter()
for ai in active_atoms:
    d = max((len(m) for m in A[ai]['poly']), default=0)
    degc[d]+=1
print(f"active atom degree dist: {dict(degc)}")
# total vars (free + gate) in the cone
allvars = set()
for ai in active_atoms:
    allvars |= L.atom_vars(A[ai]['poly'])
print(f"total distinct vars in active atoms: {len(allvars)}")
print(f"  of which free: {len(allvars & freeinp)}, gate: {len(allvars - freeinp)}")
pickle.dump({'seen':seen,'active_atoms':active_atoms,'allvars':allvars}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/cone.pkl','wb'))
