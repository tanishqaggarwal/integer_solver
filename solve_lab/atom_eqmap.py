#!/usr/bin/env python3
"""Confirm: the 3 nonzero atoms at wire=p (agentA) map to exactly the 11 failing equations.
Also check atom 37110 (forcing lock) appears in the 13 unpacking eqs and has no slack."""
import json,re
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
# Need atom->equation map. Equations reference atoms by their sub-expression. 
# Reconstruct via: an equation "uses" an atom if the atom's variable-set is a subset and ... 
# Simpler: rebuild which eqs are the 11 and confirm they involve {642,2099,7068,4432,19964,28730} (gap vars)
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
lines=open('../EQUATIONS.txt').read().split('\n')
gapvars={642,2099,7068,4432,19964,28730,23754,6947,26874,9413,17499,17325,28599}
VAR=re.compile(r'x_(\d+)')
print("=== The 11 failing eqs (wire=p) and their gap-var content ===")
for i in FAILS11:
    vs=set(int(m) for m in VAR.findall(lines[i]))
    print(f"  eq {i}: contains gap vars {sorted(vs&gapvars)}")
# atom 37110 = p - x_26064: which eqs contain x_26064 AND the constant p pattern
unp=[8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666]
print(f"\n=== forcing atom 37110 = p - x_26064; the 13 unpacking eqs all contain x_26064: ===")
cnt=sum(1 for i in unp if 'x_26064' in lines[i])
print(f"  {cnt}/13 unpacking eqs contain x_26064")
# Does x_26064 have any product-slack? check atoms containing x_26064 that are products
prodslack=[]
for ai,poly in enumerate(atoms):
    if 26064 in atom_vars(poly):
        for m in poly:
            if len(m)==2 and 26064 in m: prodslack.append(ai); break
print(f"  atoms where x_26064 appears in a PRODUCT (potential slack): {len(prodslack)} -> {prodslack[:10]}")
