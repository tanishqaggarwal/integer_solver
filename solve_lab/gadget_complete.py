import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('gadget_zeroed.json')   # already has gadget residues zeroed
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
F0=set(H.fails())
print(f"start (gadget-zeroed): {len(F0)} fails: {sorted(F0)}")
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
for n in [9629,6947,33168,950]:
    print(f"x_{n}: {'FREE' if n in H.freeinp else 'gate'} def={gdef.get(n,('',()))[0][:50]}")
print(f"x_9106={V[9106]}, x_9106/13523997={V[9106]//13523997 if V[9106]%13523997==0 else 'NOT DIV'}")
print(f"x_2239%p={V[2239]%p}, x_2239//p={V[2239]//p}")
print(f"atom17897 val = x_9106 - 13523997*x_9629 = {V[9106]-13523997*V[9629]}")
print(f"x_23754={V[23754]}, x_26874={V[26874]}, x_6947={V[6947]}")
