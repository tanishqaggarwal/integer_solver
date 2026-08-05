#!/usr/bin/env python3
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward(); base=set(H.fails())
print(f"base: {len(base)} fails; x_12553 free={12553 in H.freeinp}, val={str(vA[12553])[:10]}..")
print(f"x_19964 currently = x_12553? {vA[19964]==vA[12553]}")
# set x_12553 = CONST1 (=x_4432), so x_19964 -> CONST1, fixing G2
H.val[12553]=CONST1
H.forward(); F=set(H.fails())
print(f"\nset x_12553=CONST1: {len(F)} fails; x_19964==x_4432 now? {H.val[19964]==vA[4432]}")
print(f"  fixed: {sorted(base-F)}")
print(f"  broke: {len(F-base)} -> {sorted(F-base)[:25]}")
# G1 analog: is x_2099 = a clean free residue too? x_2099=x_37158+x_25297, x_37158=x_10878+x_22542
for lbl,v in [('x_37158',37158),('x_10878',10878),('x_22542',22542),('x_25297',25297)]:
    print(f"  {lbl} free={v in H.freeinp}, val={str(vA[v])[:8] if abs(vA[v])<10**11 else 'BIG'}")
