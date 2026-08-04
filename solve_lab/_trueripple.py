import heal_harness as H, re
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
F0=set(H.fails())
print("base fails:",len(F0))
for knob in [7068,4432,17325,9413]:
    setfree(knob, base[knob]+1)
    F1=set(H.fails())
    setfree(knob, base[knob])
    broke=sorted(F1-F0); fixed=sorted(F0-F1)
    print(f"perturb x_{knob} by +1: now {len(F1)} fails; broke {len(broke)} satisfied eqs: {broke[:30]}; fixed {fixed}")
