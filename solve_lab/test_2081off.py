import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
F0=set(H.fails())
print(f"baseline: x_2081={V[2081]}, x_4287={V[4287]}, x_15298={V[15298]}, x_21279={V[21279]}, fails={len(F0)}")
print(f"x_1222={V[1222]} x_35723={V[35723]} x_28505={V[28505]} x_32083={V[32083]}")
# set x_2081=0
V[2081]=0
H.forward()
print(f"\nafter x_2081=0: x_15298={V[15298]}, x_20434={V[20434]}, x_21279={V[21279]}")
F=set(H.fails())
print(f"fails={len(F)}")
print(f"  fixed: {sorted(F0-F)[:20]}")
print(f"  broken: {sorted(F-F0)[:30]}")
# does x_2081 feed x_15298 at all?
print(f"\nx_2081 in anc of x_15298? {2081 in H.anc.get(15298,set())}")
