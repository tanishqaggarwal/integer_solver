#!/usr/bin/env python3
import heal_harness as H
from collections import Counter
p=H.p
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F16=H.fails()
print(f"state: {len(F16)} fail: {F16}")
PIN={4432,7068,14853,31339}
# for each of the 16: free inputs present (excluding pinned), with coefficient
def coefv(i,w):
    ns={'v':H.val,'__builtins__':{}}
    r0=eval(H.eqcode[i],ns); H.val[w]+=1; r1=eval(H.eqcode[i],ns); H.val[w]-=1
    return (r1-r0)
for i in F16:
    fv=[w for w in H.eqvars[i] if w in H.freeinp and w not in PIN]
    # coefficients
    info=[]
    for w in fv:
        c=coefv(i,w)%p
        info.append((w,c))
    print(f"eq {i}: free knobs(excl pin)={[(w) for w,c in info]}")
# which free inputs appear across multiple of the 16?
allfv=Counter()
for i in F16:
    for w in H.eqvars[i]:
        if w in H.freeinp and w not in PIN: allfv[w]+=1
print("\nfree knobs by frequency across 16:", allfv.most_common(20))
