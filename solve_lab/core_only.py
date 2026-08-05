#!/usr/bin/env python3
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
F16=set([697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431])
ORIG11=set([2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125])
# 39013 core-fail set
for v in H.freeinp: H.val[v]=v013.get(v,0)
H.forward(); F013=set(H.fails())
print(f"39013: {len(F013)} fail (the core)")
# apply ONLY core knobs x_14853, x_31339
for v in H.freeinp: H.val[v]=v013.get(v,0)
H.val[14853]=vA[14853]; H.val[31339]=vA[31339]
H.forward(); Fc=set(H.fails())
print(f"\napply ONLY x_14853,x_31339 (core knobs): {len(Fc)} fail")
print(f"  core still failing? {len(Fc & F013)} of the 20 core")
print(f"  16 broken? {len(Fc & F16)} of the 16")
print(f"  orig-11 broken? {len(Fc & ORIG11)}")
print(f"  NEW fails (not core, not 16, not 11): {sorted(Fc - F013 - F16 - ORIG11)[:30]}")
print(f"  total NEW (not in 39013 core): {len(Fc - F013)}")
