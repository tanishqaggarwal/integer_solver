#!/usr/bin/env python3
"""Constrained mod-p linear heal: zero the 16 ripple eqs using a compensator pool of free inputs,
while keeping all touched equations satisfied. Pin x_4432,x_7068 (protect the 11) + core knobs."""
import heal_harness as H
import time
p=H.p
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F16=set(H.fails())
print(f"start: {len(F16)} fail")
PIN={4432,7068,14853,31339}
# Pool: all free inputs in the 16 (excl pins) + also include descendants' free ancestors
pool=set()
for i in F16:
    for w in H.eqvars[i]:
        if w in H.freeinp and w not in PIN: pool.add(w)
pool=sorted(pool)
print(f"compensator pool size: {len(pool)}")
# descendants of pool -> equations touched
t0=time.time()
desc=set(pool)
for t in H.order:
    if H.anc[t]&set(pool): desc.add(t)
Etouch=[i for i in range(len(H.eqcode)) if H.eqvars[i]&desc]
print(f"descendant gate outputs: {len(desc)}, equations touched: {len(Etouch)} ({time.time()-t0:.1f}s)")
# baseline residuals of Etouch
ns={'v':H.val,'__builtins__':{}}
r0={i:eval(H.eqcode[i],ns) for i in Etouch}
# Jacobian: bump each pool var by 1, forward-reconstruct, measure delta resid on Etouch
J={}  # (row_i, col_w) -> coef
base=H.val[:]
for w in pool:
    H.val[w]+=1
    H.forward()
    ns={'v':H.val,'__builtins__':{}}
    for i in Etouch:
        d=(eval(H.eqcode[i],ns)-r0[i])
        if d%p!=0: J[(i,w)]=d%p
    H.val[:]=base[:]  # restore
    H.forward()
print(f"Jacobian nnz: {len(J)} ({time.time()-t0:.1f}s)")
# Build target: 16 -> -r0, others -> 0
import pickle
pickle.dump({'pool':pool,'Etouch':Etouch,'r0':{i:r0[i]%p for i in Etouch},'J':J,'F16':sorted(F16)}, open('heal16.pkl','wb'))
print("saved heal16.pkl")
