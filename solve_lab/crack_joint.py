#!/usr/bin/env python3
import heal_harness as H
from math import gcd
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward(); base=set(H.fails())
print(f"base: {len(base)} fails: {sorted(base)}")
# gaps: x_19964 should become x_4432; x_2099 should become x_7068 (so G1,G2 hold with slacks 0)
# x_19964 = x_20492 + x_21279*x_8731 ; want = x_4432 => x_21279*x_8731 = x_4432 - x_20492
# x_2099  = x_37158 + x_21279*x_9118 ; want = x_7068 => x_21279*x_9118 = x_7068 - x_37158
g2 = vA[4432]-vA[20492]
g1 = vA[7068]-vA[37158]
print(f"g2 (x_21279*x_8731) = {g2}")
print(f"g1 (x_21279*x_9118) = {g1}")
g = gcd(g2,g1)
print(f"gcd(g2,g1) = {g}  (x_21279=gcd, x_8731=g2/gcd, x_9118=g1/gcd)")
# x_21279 = x_4287 (since x_9062=x_4287, x_20434=1). Set x_4287=g, x_8731=g2//g, x_9118=g1//g
H.val[4287]=g; H.val[8731]=g2//g; H.val[9118]=g1//g
H.forward()
F=set(H.fails())
print(f"\nx_4287={g<10**8 and g or 'big'}, x_8731={g2//g if abs(g2//g)<10**8 else 'big'}, x_9118={g1//g if abs(g1//g)<10**8 else 'big'}")
print(f"after joint fix: {len(F)} fails")
print(f"  x_19964==x_4432? {H.val[19964]==vA[4432]}; x_2099==x_7068? {H.val[2099]==vA[7068]}")
print(f"  fixed: {sorted(base-F)}")
print(f"  broke: {len(F-base)} -> {sorted(F-base)[:25]}")
