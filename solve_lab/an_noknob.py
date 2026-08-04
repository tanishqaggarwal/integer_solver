#!/usr/bin/env python3
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
# Which changed frees are ancestors of core quantities S=x_35389, T=x_6671, loads?
core_targets=[35389,6671,11150,25739,37758,15298,29804,35605,4007,29322,3558,24908,14853,16742,12186,33469,27713,29356,1326]
core_anc=set()
for ct in core_targets: core_anc |= H.anc.get(ct,{ct})
print("changed frees that are core ancestors:", [v for v in changed_free if v in core_anc])
print("changed frees NOT core ancestors:", [v for v in changed_free if v not in core_anc])

# Test: apply all changed frees EXCEPT 4432,7068
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F=H.fails()
print(f"\napply 11 knobs (skip 4432,7068), forward: {len(F)} fail")
# is core among fails? core squares - check the 20. Use known core detection: eqs that were failing at 39013
F013=set()
for v in H.freeinp: H.val[v]=v013.get(v,0)
H.forward(); F013=set(H.fails())
print(f"39013 fails (the core): {len(F013)}: {sorted(F013)}")
# recompute the test state
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward(); Ft=set(H.fails())
print(f"\ntest-state fails ({len(Ft)}): {sorted(Ft)}")
print(f"  core still failing (of the 39013 core): {sorted(Ft&F013)}")
print(f"  new (ripple) fails not in core: {sorted(Ft-F013)}")
print(f"  original-11 among fails: {sorted(Ft&{2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125})}")
