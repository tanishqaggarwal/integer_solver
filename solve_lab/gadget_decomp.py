#!/usr/bin/env python3
import heal_harness as H
import json
p=H.p
d=json.load(open('fc_partial.json'))
fp={int(k[2:]) if k.startswith('x_') else int(k):int(v) for k,v in d.items()}
for v in H.freeinp: H.val[v]=fp.get(v,0)
H.forward()
import json as J
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=J.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
def g(t): return gdef.get(t,('FREE',()))[0]
# structure
for t in [27019,37530,37944,21671,25490,27177,4306,7181,7010,3349,17925,34310]:
    print(f"x_{t} = {g(t)[:40]}  val={'0' if H.val[t]==0 else 'nz'} free={t in H.freeinp}")
# Is x_27019 linear in x_14865? evaluate at 0,1,2
def ev27019(x14865):
    old=H.val[14865]; H.val[14865]=x14865; H.forward(); r=H.val[27019]; H.val[14865]=old; H.forward(); return r
a=ev27019(0); b=ev27019(1); c=ev27019(2)
print(f"\nx_27019 at x_14865=0,1,2: linear? 2nd-diff={c-2*b+a} (0=linear)")
print(f"  slope={b-a}, so x_27019=0 at x_14865 = {(-(a))//(b-a) if (b-a)!=0 and a%(b-a)==0 else 'noninteger/'+str(a)+'/'+str(b-a)}")
