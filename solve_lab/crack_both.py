#!/usr/bin/env python3
import heal_harness as H
import pickle
p=H.p
CK=pickle.load(open('checked.pkl','rb')); checked=CK['checked']
vA=H.loadd('best_agentA_39022.json')
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward(); base=set(H.fails())
print(f"base: {len(base)} fails")
print(f"x_12553 checked={12553 in checked}, x_6418 checked={6418 in checked}")
H.val[12553]=CONST1; H.val[6418]=CONST2
H.forward(); F=set(H.fails())
print(f"set x_12553=CONST1, x_6418=CONST2: {len(F)} fails")
print(f"  x_19964==x_4432? {H.val[19964]==vA[4432]}; x_2099==x_7068? {H.val[2099]==vA[7068]}")
print(f"  fixed (of orig 11): {sorted(base&(base-F))} -> {sorted(base-F)}")
print(f"  broke: {len(F-base)} -> {sorted(F-base)}")
# Are the broken ones checks on x_12553/x_6418 or residue-consumers? show what they need
print(f"\nremaining fails: {sorted(F)}")
