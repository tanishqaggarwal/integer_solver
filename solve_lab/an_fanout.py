#!/usr/bin/env python3
import heal_harness as H
from collections import defaultdict
# fan-out: number of equations each free input touches (via descendants)
# precompute eq -> free ancestors
eq_free=[]
for i in range(len(H.eqcode)):
    s=set()
    for w in H.eqvars[i]:
        if w in H.freeinp: s.add(w)
        s|=H.anc.get(w,set())
    eq_free.append(s)
free_eqs=defaultdict(int)
for i in range(len(H.eqcode)):
    for w in eq_free[i]: free_eqs[w]+=1
# core difference variables and their knobs
for label,v in [('x_14853',14853),('x_12186',12186),('x_16742',16742),('x_24908',24908),('x_31339',31339),
                ('x_1326',1326),('x_29322',29322),('x_3558',3558)]:
    isfree = v in H.freeinp
    print(f"{label}: free={isfree}, touches {free_eqs.get(v,'(gate-out)') if isfree else 'N/A'} eqs" + ("" if isfree else f"; free-ancestors={len(H.anc.get(v,set()))}"))
# Also: what free inputs feed x_24908 (for x_3558)?
print("\nx_24908 free ancestors:", sorted(H.anc.get(24908,set()))[:20], "...total", len(H.anc.get(24908,set())))
# fan-out of x_24908's ancestors
for w in sorted(H.anc.get(24908,set())):
    print(f"  x_{w}: touches {free_eqs.get(w,0)} eqs")
