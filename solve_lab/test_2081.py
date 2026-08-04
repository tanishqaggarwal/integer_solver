#!/usr/bin/env python3
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for cfg in [{}, {2081:0}, {2081:0, 8599:1, 25956:1}]:
    for v in H.freeinp: H.val[v]=vA.get(v,0)
    for k,val in cfg.items(): H.val[k]=val
    H.forward()
    F=H.fails()
    print(f"cfg={cfg}: x_15298={H.val[15298]}, x_5814={H.val[5814]}, fails={len(F)}")
