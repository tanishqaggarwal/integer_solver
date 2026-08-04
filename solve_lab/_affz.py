import heal_harness as H, json, re
from collections import defaultdict
p=H.p
VAR=re.compile(r'x_(\d+)')
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
# union of free vars in failing eqs
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
relfree=set()
for i in F:
    lhs=lines[i].rsplit('=',1)[0]
    for m in VAR.findall(lhs):
        v=int(m)
        if v in H.freeinp: relfree.add(v)
relfree=sorted(relfree)
print("relevant free vars (%d):"%len(relfree),relfree)
# desc_of for incremental forward
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def resid():
    return [eval(H.eqcode[i],ns) for i in F]
# test affinity in Z: second difference for each free var
r0=resid()
affine_in={}
for v in relfree:
    setfree(v, base[v]+1); r1=resid()
    setfree(v, base[v]+2); r2=resid()
    setfree(v, base[v])   # restore
    # second diff per eq
    d1=[r1[k]-r0[k] for k in range(len(F))]
    d2=[r2[k]-r1[k] for k in range(len(F))]
    isaff = all(d2[k]==d1[k] for k in range(len(F)))
    affine_in[v]=(isaff, d1)
    tag="AFFINE" if isaff else "NONLINEAR"
    nz=[F[k] for k in range(len(F)) if d1[k]!=0 or d2[k]!=0]
    print(f"  x_{v}: {tag}  affects eqs {nz}")
