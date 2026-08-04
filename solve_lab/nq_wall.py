import heal_harness as H
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
print(f"x_2239 % p = {V[2239]%p}")
print(f"x_2239 is {'ZERO mod p' if V[2239]%p==0 else 'NONZERO mod p (potential wall)'}")
# x_2239 free ancestors
anc2239=H.anc.get(2239,set())&H.freeinp
print(f"x_2239 free ancestors ({len(anc2239)}): {sorted(anc2239)}")
# find x_2239 gate def
gdef={}
import json
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
print(f"x_2239 def: {gdef.get(2239,('?',()))[0][:200]}")
print(f"x_23754 def: {gdef.get(23754,('?',()))[0][:200]}")
print(f"x_26874 def: {gdef.get(26874,('?',()))[0][:200]}")
# is x_26874 a wire member (=p)?
print(f"x_26874 % p = {V[26874]%p}, ==p? {V[26874]==p}")
# check: can moving x_2239's ancestors change x_2239 mod p?
# perturb each free ancestor by 1, see effect on x_2239 mod p
ns={'v':V,'__builtins__':{}}
from collections import defaultdict
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
base=V[2239]%p
print(f"\nsensitivity of x_2239 mod p to ancestors:")
for w in sorted(anc2239):
    old=V[w]; V[w]=old+1
    for k in desc_of[w]:
        if H.order[k]==2239 or 2239 in [H.order[k]]: pass
    for k in desc_of[w]: V[H.order[k]]=eval(H.gcode[k],ns)
    delta=(V[2239]%p - base)%p
    V[w]=old
    for k in desc_of[w]: V[H.order[k]]=eval(H.gcode[k],ns)
    print(f"  x_{w}: d(x_2239 mod p)/dx = {delta if delta<10**6 or delta>p-10**6 else str(delta)[:20]+'...'}")
