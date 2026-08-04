#!/usr/bin/env python3
import heal_harness as H
import json
p=H.p
d=json.load(open('fc_partial.json'))
fp={int(k[2:]) if k.startswith('x_') else int(k):int(v) for k,v in d.items()}
for v in H.freeinp: H.val[v]=fp.get(v,0)
H.forward()
F0=set(H.fails())
print(f"base: {len(F0)} fails")
# THE FIX: x_14865=x_12553, x_31861=x_6418, x_33168=0
H.val[14865]=H.val[12553]
H.val[31861]=H.val[6418]
H.val[33168]=0
H.forward()
F1=set(H.fails())
print(f"after x_14865=x_12553, x_31861=x_6418, x_33168=0: {len(F1)} fails")
print(f"  x_27177={H.val[27177]}, x_4306={H.val[4306]}, x_31731={H.val[31731]}, x_9106%13523997={H.val[9106]%13523997}")
print(f"  fixed: {sorted(F0-F1)}")
print(f"  broke: {sorted(F1-F0)}")
if len(F1)==0:
    json.dump({f'x_{i}':H.val[i] for i in range(H.NVARS)},open('SOLVED_full.json','w'))
    print("*** SOLVED 39033/39033 - saved SOLVED_full.json ***")
elif len(F1)<len(F0):
    json.dump({f'x_{i}':H.val[i] for i in range(H.NVARS)},open('gadget_progress.json','w'))
    print(f"progress: {39033-len(F1)}/39033, saved gadget_progress.json")
