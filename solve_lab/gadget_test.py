#!/usr/bin/env python3
import heal_harness as H
import json
p=H.p
d=json.load(open('fc_partial.json'))
fp={int(k[2:]) if k.startswith('x_') else int(k):int(v) for k,v in d.items()}
for v in H.freeinp: H.val[v]=fp.get(v,0)
H.forward()
print(f"base: {len(H.fails())} fails; x_27177%p={H.val[27177]%p!=0}, x_4306%p={H.val[4306]%p!=0}")
# sensitivity of x_27177, x_4306, x_9106 to x_14865, x_31861 (finite diff)
def probe(w):
    b27,b43,b91=H.val[27177],H.val[4306],H.val[9106]
    H.val[w]+=1; H.forward()
    d27=(H.val[27177]-b27)%p; d43=(H.val[4306]-b43)%p; d91=(H.val[9106]-b91)%p
    H.val[w]-=1; H.forward()
    return d27,d43,d91
for w in [14865,31861,6418,12553,9118,8731]:
    d27,d43,d91=probe(w)
    print(f"  d/dx_{w}: x_27177{'+' if d27 else '0'} x_4306{'+' if d43 else '0'} x_9106{'+' if d91 else '0'}  (27177 lin={d27}, 9106 lin={d91%13523997})")
