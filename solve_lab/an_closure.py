#!/usr/bin/env python3
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F16=set(H.fails())
PIN={4432,7068,14853,31339}
# closure: alternate pool(free) <-> touched(eqs)
# eq -> free ancestors ; free -> eqs it touches (via descendants)
# precompute free-ancestors per eq
eq_free=[set(w for w in H.eqvars[i] if w in H.freeinp) | set().union(*[H.anc[w] for w in H.eqvars[i] if w in H.anc]) if H.eqvars[i] else set() for i in range(len(H.eqcode))]
# simpler: eq_free[i] = union of anc of each var in eq
eq_free=[]
for i in range(len(H.eqcode)):
    s=set()
    for w in H.eqvars[i]:
        if w in H.freeinp: s.add(w)
        s|=H.anc.get(w,set())
    eq_free.append(s)
# free -> eqs
from collections import defaultdict
free_eqs=defaultdict(set)
for i in range(len(H.eqcode)):
    for w in eq_free[i]:
        free_eqs[w].add(i)
pool=set(w for i in F16 for w in eq_free[i])-PIN
eqs=set(F16)
for it in range(40):
    newe=set()
    for w in pool: newe|=free_eqs[w]
    newf=set()
    for i in newe: newf|=eq_free[i]
    newf-=PIN
    if newe<=eqs and newf<=pool:
        break
    eqs|=newe; pool|=newf
print(f"closure after {it} iters: pool(free)={len(pool)}, eqs={len(eqs)}")
print(f"(total free inputs={len(H.freeinp)}, total eqs={len(H.eqcode)})")
