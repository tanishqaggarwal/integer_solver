#!/usr/bin/env python3
"""Classify atoms; report equation-level statistics."""
import pickle, re, os, sys, collections, ast

HERE = os.path.dirname(os.path.abspath(__file__))
D = pickle.load(open(os.path.join(HERE, 'atoms.pkl'), 'rb'))
atom_src = D['atom_src']; atom_vars = D['atom_vars']; eq_terms = D['eq_terms']

print("n_atoms", len(atom_src), "n_eqs", len(eq_terms))

# equations by number of atoms
c = collections.Counter(len(t) for t in eq_terms)
print("eq atom-count histogram (first 12):", sorted(c.items())[:12])

single = [i for i, t in enumerate(eq_terms) if len(t) == 1]
print("single-atom equations:", len(single))
forced_atoms = set(t[0][1] for t in (eq_terms[i] for i in single))
print("atoms forced to zero by single-atom eqs:", len(forced_atoms))

# classify atom shapes
pat_def = re.compile(r'^X(\d+) - ')
shapes = collections.Counter()
for s in atom_src:
    t = re.sub(r'X\d+', 'V', s)
    t = re.sub(r'\d+', 'N', t)
    shapes[t] += 1
print("\ntop 40 atom shapes:")
for sh, k in shapes.most_common(40):
    print(f"  {k:6d}  {sh}")
print("distinct shapes:", len(shapes))
