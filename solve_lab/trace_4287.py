#!/usr/bin/env python3
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
gdef={}; consumers=defaultdict(list)
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gdef[d['t']]=(d['rhs'],tuple(d['vids']))
        for u in d['vids']: consumers[u].append(d['t'])
# x_4287 direct consumers (gates using it)
print(f"x_4287 direct gate consumers: {consumers[4287]}")
for c in consumers[4287]:
    print(f"   x_{c} = {gdef[c][0][:50]}")
# x_4287 in atoms directly?
inat=[ai for ai,pp in enumerate(atoms) if 4287 in atom_vars(pp)]
print(f"x_4287 in {len(inat)} atoms directly")
# x_8731, x_9118 consumers/atoms (bits)
for b in [8731,9118]:
    inatb=[ai for ai,pp in enumerate(atoms) if b in atom_vars(pp)]
    # is there a boolean-enforcement atom x_b*(x_b-1)? i.e., atom with monomial (b,b)
    boolatom=[ai for ai in inatb if (b,b) in atoms[ai]]
    print(f"x_{b}: in {len(inatb)} atoms, consumers={consumers[b][:8]}, has x*(x-1) boolean atom: {boolatom}")
