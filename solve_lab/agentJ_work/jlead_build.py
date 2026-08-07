#!/usr/bin/env python3
"""Rebuild jlead.pkl (syntactic definer candidate per atom) from jmodel2/jpoly."""
import pickle, re, os
HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
atoms, dc = M['atoms'], P['defcands']
pat = re.compile(r'^\(+x_(\d+)\)+-')
lead = [None] * len(atoms)
for i, a in enumerate(atoms):
    m = pat.match(a)
    if m:
        v = int(m.group(1))
        if v in dc[i]:
            lead[i] = v
pickle.dump(lead, open(os.path.join(HERE, 'jlead.pkl'), 'wb'))
print("lead definers:", sum(1 for v in lead if v is not None))
