#!/usr/bin/env python3
"""Check whether the 9 monster free inputs appear in equations OTHER than the 20 core.
If private, we can freely set them to zero M1,M2,M3."""
import json, re, sys
from propagate import NVARS
p=2**256-2**32-977
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
VAR=re.compile(r'x_(\d+)')
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
monster_frees=[21889,32405,30317, 25156,3023,5146, 2287,15663,2936]
print("free input -> equations containing it (C=core, L=linear-satisfied):")
for v in monster_frees:
    occ=[i for i in range(len(lines)) if v in eqvars[i]]
    core_occ=[i for i in occ if i in CORE]
    lin_occ=[i for i in occ if i not in CORE]
    print(f"  x_{v}: total {len(occ)} eqs | core={core_occ} | linear={lin_occ[:15]}{'...' if len(lin_occ)>15 else ''}")
