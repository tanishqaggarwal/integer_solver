#!/usr/bin/env python3
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
base=set(H.fails())
print(f"agentA base: {len(base)} fails; free? x_4287={4287 in H.freeinp}, x_8731={8731 in H.freeinp}")
# target: x_19892 = x_4432 - x_20492 so that x_19964 = x_20492 + x_19892 = x_4432
gap = vA[4432]-vA[20492]
print(f"gap (need x_19892) = {gap}")
# x_21279 = x_9062*x_20434 = x_4287*x_20434. x_20434=1 => x_21279=x_4287. x_19892=x_21279*x_8731=x_4287*x_8731.
# set x_8731=1, x_4287=gap
H.val[4287]=gap; H.val[8731]=1
H.forward()
F=set(H.fails())
print(f"\nset x_4287=gap, x_8731=1, forward: {len(F)} fails")
print(f"  x_19964 now = x_4432? {H.val[19964]==vA[4432]}  (x_19964={str(H.val[19964])[:12]}..)")
print(f"  x_21279 now = {str(H.val[21279])[:12]}.., x_19892={str(H.val[19892])[:12]}..")
print(f"  fixed: {sorted(base-F)}")
print(f"  broke: {len(F-base)} -> {sorted(F-base)[:30]}")
