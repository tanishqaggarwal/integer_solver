#!/usr/bin/env python3
"""Close the final 12. They need x_33462=CONST1, x_22152=CONST2 (loads corrupted by cascade).
Diagnose: (A) set them to CONST; (B) turn off activator x_24601; inspect x_24601 role."""
import json
from agentA_harness import (p, load_solution, forward, backward_cone, freeinp, eqcode,
                            eqvars, lines, NEQ, NVARS, gates, definer, order)
CONST1 = 97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2 = 126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
CORE = {2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}
TWELVE = [4833,4944,5348,9344,10406,11574,12321,19708,20927,21972,27514,38014]

v0 = load_solution('best_agentA_39021.json'); forward(v0)
def count(v):
    ns={'__builtins__':{},'v':v}; return set(i for i in range(NEQ) if eval(eqcode[i],ns)!=0)
F0 = count(v0)
print(f"39021 baseline: {NEQ-len(F0)}/{NEQ} ({len(F0)} fail); core fails={len(F0&CORE)}")
print(f"current x_22152={v0[22152]}  (==CONST2? {v0[22152]==CONST2})")
print(f"current x_33462={v0[33462]}  (==CONST1? {v0[33462]==CONST1})")
print(f"x_22152 free? {22152 in freeinp}  x_33462 free? {33462 in freeinp}  x_24601 free? {24601 in freeinp}")
print(f"x_24601 value={v0[24601]}  in #eqs={sum(1 for i in range(NEQ) if 24601 in eqvars[i])}")

# does x_24601 feed x_15298?
_,c15298 = backward_cone(15298)
print(f"x_24601 feeds x_15298? {24601 in c15298}   x_2081 feeds x_15298? {2081 in c15298}")
print(f"x_15298 current = {v0[15298]}")

# Option A: set loads to CONST
vA = v0[:]; vA[22152]=CONST2; vA[33462]=CONST1; forward(vA)
FA = count(vA)
print(f"\n[A] set x_22152=CONST2,x_33462=CONST1: {NEQ-len(FA)}/{NEQ} ({len(FA)} fail)")
print(f"    fixed(of 12): {sorted(set(TWELVE)-FA)}  still-fail-12: {sorted(set(TWELVE)&FA)}")
print(f"    broke(new): {sorted(FA-F0)[:30]}  core fails={len(FA&CORE)}")

# Option B: turn off x_24601
vB = v0[:]; vB[24601]=0; forward(vB)
FB = count(vB)
print(f"\n[B] set x_24601=0: {NEQ-len(FB)}/{NEQ} ({len(FB)} fail); x_15298={vB[15298]}")
print(f"    fixed(of 12): {sorted(set(TWELVE)-FB)}  still-fail-12: {sorted(set(TWELVE)&FB)}")
print(f"    broke(new): {sorted(FB-F0)[:30]}  core fails={len(FB&CORE)}")

# what does x_24601 multiply in the 12? and elsewhere
import re
cnt=0
for i in range(NEQ):
    if 24601 in eqvars[i]: cnt+=1
print(f"\nx_24601 appears in {cnt} eqs. Sample contexts:")
seen=set()
for i in range(NEQ):
    if 24601 in eqvars[i]:
        for m in re.finditer(r'.{4}x_24601.{22}', lines[i]):
            s=m.group(0)
            if s not in seen: seen.add(s); print(f"   {s}")
        if len(seen)>=8: break
