import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
V=H.val
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
for n in [7075,19247,642,28730,17325,9413,2099,19964,6418,12553]:
    print(f"x_{n}: {'FREE' if n in H.freeinp else 'gate'} def={gdef.get(n,('',()))[0][:45]}")
V[2081]=0
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
# zero x_7075, x_19247 via free ancestors
print(f"\nx_7075 free-anc: {sorted(H.anc.get(7075,set())&H.freeinp)}")
print(f"x_19247 free-anc: {sorted(H.anc.get(19247,set())&H.freeinp)}")
print(f"x_2099={V[2099]}, x_9118={V[9118]}, equal={V[2099]==V[9118]}")
print(f"x_19964={V[19964]}, x_8731={V[8731]}, equal={V[19964]==V[8731]}")
F0=set(H.fails())
print(f"\nx_2081=0 baseline fails: {len(F0)}")
