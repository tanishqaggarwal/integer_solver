import json
import heal_harness as H
from collections import defaultdict
p=H.p
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
# footprint of each free input
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
eq_frees=[]
for i in range(len(H.lines)):
    s=set()
    for v in H.eqvars[i]: s|=H.anc.get(v, {v} if v in H.freeinp else set())
    eq_frees.append(s & H.freeinp)
foot=defaultdict(set)
for i,fs in enumerate(eq_frees):
    for f in fs: foot[f].add(i)
for n in [19083,1308,9254,29967,16742,14853,12186,24908]:
    fr='FREE' if n in H.freeinp else 'gate'
    df=gdef.get(n,('',()))[0][:40]
    ft=len(foot.get(n,set())) if n in H.freeinp else '-'
    print(f"x_{n}: {fr} def={df} footprint={ft} eqs")
# key: to heal atom27902 (12846437*(x_14853-x_1308)-x_29967) after x_14853 changed,
# keep x_14853-x_1308 constant by moving x_1308, OR adjust x_29967.
# Which is free and low-footprint?
print("\nInterpretation:")
print(" atom27902: x_29967 = 12846437*(x_14853 - x_1308)")
print(" atom25170: x_9254  = 6788513*(x_16742 - x_19083)")
