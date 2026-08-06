#!/usr/bin/env python3
"""Test whether the 16 ripple + other touched equations are LINEAR in the compensator free inputs.
Bump a free by 1 and by a large L; if residual change scales exactly by L, it's linear in that free."""
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
ns={'v':H.val,'__builtins__':{}}
def R(idxs): return {i:eval(H.eqcode[i],ns) for i in idxs}
base=R(F16)
# test frees: compensators + core knobs
tests=[24548,2964,9280,27711,15324,14853,31339,11080,4432,7068]
L=10**40
import random
print("free : linear-in-16? (all 16 residuals scale exactly with bump)")
for w in tests:
    # bump by 1
    H.val[w]+=1; H.forward(); r1=R(F16)
    H.val[w]+=(L-1); H.forward(); rL=R(F16)
    H.val[w]-=L; H.forward()
    # linear iff (rL[i]-base[i]) == L*(r1[i]-base[i]) for all i
    lin=all((rL[i]-base[i])==L*(r1[i]-base[i]) for i in F16)
    nz=sum(1 for i in F16 if r1[i]!=base[i])
    print(f"  x_{w}: affects {nz}/16, LINEAR={lin}")
