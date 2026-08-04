#!/usr/bin/env python3
import heal_harness as H
import json
from collections import defaultdict
p=H.p
# load fc_partial free inputs
d=json.load(open('fc_partial.json'))
fp={int(k[2:]) if k.startswith('x_') else int(k):int(v) for k,v in d.items()}
for v in H.freeinp: H.val[v]=fp.get(v,0)
H.forward()
F0=set(H.fails())
print(f"fc_partial reconstructed: {len(F0)} fails")
# The 6 gadget targets: x_27177=0, x_4306=0 (=> x_31731=0, x_9106=0 div), x_33168=0
print(f"current: x_27177={H.val[27177]!=0 and 'nonzero' or 0}, x_4306={H.val[4306]!=0 and 'nonzero' or 0}, x_31731={H.val[31731]!=0 and 'nz' or 0}, x_9106%13523997={H.val[9106]%13523997}, x_33168={H.val[33168]}")
# x_27177 = x_21671 - x_37530 = x_7181*x_7010 - x_27019^2 ; x_4306 = x_25490 - x_37944 = x_3349*x_17925 - x_27019*x_34310
# free inputs controlling these
for t in [27177,4306,21671,37530,25490,37944,27019]:
    fa=sorted(H.anc.get(t,set())&H.freeinp)
    print(f"  x_{t}: free-ancestors={fa[:10]}{'...' if len(fa)>10 else ''} (n={len(fa)})")
