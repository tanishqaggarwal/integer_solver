#!/usr/bin/env python3
"""Pass 3: atoms / equations supported on BITS + exactly one extra variable."""
import pickle, json, os, sys
from collections import Counter, defaultdict
from scan import load, support, degree, fmt

HERE = os.path.dirname(os.path.abspath(__file__))

D, chain, BITS = load()
atoms = D['atoms']; eq_terms = D['eq_terms']; eq_poly = D['eq_poly']
eq_outer = D['eq_outer']
idx = {v: i for i, v in enumerate(chain)}

atom_sup = [support(a) for a in atoms]
atom_eqs = defaultdict(list)
for i, terms in enumerate(eq_terms):
    for c, aid in terms:
        atom_eqs[aid].append(i)

onestep_atoms = [i for i, s in enumerate(atom_sup)
                 if (s & BITS) and len(s - BITS) == 1]

print("=== one-step ATOMS: shape histogram (nbits, deg, nterms) ===")
sh = defaultdict(list)
for aid in onestep_atoms:
    a = atoms[aid]; s = atom_sup[aid]
    sh[(len(s & BITS), degree(a), len(a))].append(aid)
for k in sorted(sh):
    print(f"  nbits={k[0]} deg={k[1]} nterms={k[2]} -> {len(sh[k])}")

print("\n=== one-step ATOMS: representative per shape ===")
for k in sorted(sh):
    aid = sh[k][0]
    print(f"  shape {k} n={len(sh[k])}  atom#{aid} eqs={atom_eqs[aid][:5]}")
    print(f"     {fmt(atoms[aid])[:400]}")

print("\n=== one-step EQUATIONS ===")
for i, p in enumerate(eq_poly):
    s = support(p)
    if (s & BITS) and len(s - BITS) == 1:
        extra = sorted(s - BITS)[0]
        print(f"eq#{i} nbits={len(s&BITS)} extra=x_{extra} deg={degree(p)} "
              f"nterms={len(p)} outer={eq_outer[i]} natoms={len(eq_terms[i])}")
        print("     ", fmt(p)[:900])
