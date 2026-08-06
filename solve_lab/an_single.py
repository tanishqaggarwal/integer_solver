#!/usr/bin/env python3
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward(); base=set(H.fails())
print(f"agentA base: {len(base)} fail")
# try reverting each of x_4432, x_7068 individually, and combos of obstruction knobs
tests={
 'revert x_4432 only':[(4432,v013[4432])],
 'revert x_7068 only':[(7068,v013[7068])],
 'revert both':[(4432,v013[4432]),(7068,v013[7068])],
}
for name,changes in tests.items():
    for v in H.freeinp: H.val[v]=vA.get(v,0)
    for idx,val in changes: H.val[idx]=val
    H.forward(); F=set(H.fails())
    print(f"{name}: {len(F)} fail  (fixed {len(base-F)}, broke {len(F-base)})")
# Also: activate the rare-partner slacks x_9413, x_17325 to see effect
for v in H.freeinp: H.val[v]=vA.get(v,0)
# set x_9413 so x_28730 = x_4432 - x_19964 ... but x_28730=p*x_9413, need multiple of p. skip.
