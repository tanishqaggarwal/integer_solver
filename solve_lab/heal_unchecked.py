#!/usr/bin/env python3
"""Which frees affect the 16 ripple, and are they checked/unchecked? Can unchecked frees heal them?"""
import heal_harness as H
import pickle
from collections import defaultdict
p=H.p
CK=pickle.load(open('checked.pkl','rb')); checked=CK['checked']
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
ns={'v':H.val,'__builtins__':{}}
# free-ancestors of the 16
eq_anc=set()
for i in F16:
    for w in H.eqvars[i]:
        if w in H.freeinp: eq_anc.add(w)
        eq_anc|=H.anc.get(w,set())
eq_anc&=H.freeinp
chk=[w for w in eq_anc if w in checked]
unchk=[w for w in eq_anc if w not in checked]
print(f"16 ripple free-ancestors: {len(eq_anc)} total; checked={len(chk)}, unchecked={len(unchk)}")
# leverage: which unchecked frees actually change the 16 residuals?
base={i:eval(H.eqcode[i],ns) for i in F16}
lever_unchk=[]
import time; t0=time.time()
for w in unchk:
    H.val[w]+=1
    # incremental: recompute descendants
    for k,t in enumerate(H.order):
        if w in H.anc.get(t,()): H.val[t]=eval(H.gcode[k],ns)
    nz=sum(1 for i in F16 if eval(H.eqcode[i],ns)!=base[i])
    H.val[w]-=1
    for k,t in enumerate(H.order):
        if w in H.anc.get(t,()): H.val[t]=eval(H.gcode[k],ns)
    if nz>0: lever_unchk.append((w,nz))
    if time.time()-t0>60: print("(timeout partial)"); break
print(f"unchecked frees with LEVERAGE on the 16: {len(lever_unchk)}")
print("  top:", sorted(lever_unchk,key=lambda t:-t[1])[:20])
