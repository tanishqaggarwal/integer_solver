#!/usr/bin/env python3
import heal_harness as H
from collections import defaultdict
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
# eq -> free ancestors
eq_free=[]
for i in range(len(H.eqcode)):
    s=set()
    for wv in H.eqvars[i]:
        if wv in H.freeinp: s.add(wv)
        s|=H.anc.get(wv,set())
    eq_free.append(s)
# The 16's own free-ancestor closure (honest full closure), tracking eqs too
free_eqs=defaultdict(set)
for i in range(len(H.eqcode)):
    for w in eq_free[i]: free_eqs[w].add(i)
pool=set().union(*[eq_free[i] for i in F16])-PIN
eqs=set(F16)
sizes=[]
for it in range(60):
    ne=set()
    for w in pool: ne|=free_eqs[w]
    nf=set()
    for i in ne: nf|=eq_free[i]
    nf-=PIN
    if ne<=eqs and nf<=pool: break
    eqs|=ne; pool|=nf
    sizes.append((len(pool),len(eqs)))
print(f"HONEST closure of the 16: {len(pool)} free inputs, {len(eqs)} eqs (iters {len(sizes)})")
print(f"growth: {sizes[:6]} ... {sizes[-3:] if len(sizes)>3 else ''}")
# So heal_grow's '91 stopped' was because it only added ancestors of first-80 inconsistent rows.
# The TRUE compensator closure is large. Confirm the heal must span it.
print(f"\n=> heal_grow's plateau at 91 was premature (only added inconsistent-row ancestors).")
print(f"   True closure requires ~{len(pool)} free inputs / {len(eqs)} eqs -> whole-system heal.")
