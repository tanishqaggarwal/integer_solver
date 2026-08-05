import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
base[4287]=1; base[2081]=1; base[9118]=base[7068]; base[8731]=base[4432]
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F=H.fails()
print(f"branch B fails: {len(F)}")
Kset=sorted(set(v for i in F for v in H.eqvars[i] if v in H.freeinp))
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
touch={v:sorted(set(desc_of[v])) for v in Kset}
def setf(v,x):
    H.val[v]=x
    for k in touch[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def residZ():
    return {i:eval(H.eqcode[i],ns) for i in F}
r0=residZ()
nonlin_var=set(); nonlin_eq=set()
for v in Kset:
    o=base[v]
    setf(v,o+1); r1=residZ(); setf(v,o+2); r2=residZ(); setf(v,o)
    for i in F:
        if (r2[i]-r1[i])!=(r1[i]-r0[i]): nonlin_var.add(v); nonlin_eq.add(i)
print(f"nonlinear vars: {len(nonlin_var)}/{len(Kset)}")
print(f"nonlinear eqs: {len(nonlin_eq)}/{len(F)}: {sorted(nonlin_eq)}")
linear_eq=[i for i in F if i not in nonlin_eq]
print(f"LINEAR eqs: {len(linear_eq)}: {linear_eq}")
# For linear eqs, are they affine in ALL Kset (jointly)? build model & test
import random
random.seed(9)
lin_relevant=sorted(set(v for i in linear_eq for v in H.eqvars[i] if v in H.freeinp))
print(f"free vars in linear eqs: {len(lin_relevant)}")
