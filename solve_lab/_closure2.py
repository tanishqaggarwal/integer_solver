import heal_harness as H, re, random, json, time
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
def setmany(assign, touch):
    for v,x in assign.items(): H.val[v]=x
    for k in touch: H.val[H.order[k]]=eval(H.gcode[k],ns)
# var -> equations
var2eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var2eqs[v].append(i)
F0=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
# equations potentially affected by moving free var v: those containing v or a descendant of v
def affected_eqs(v):
    S=set(var2eqs.get(v,[]))
    for k in desc_of[v]: S.update(var2eqs.get(H.order[k],[]))
    return S
# ripple of x_7068,x_4432 (only evaluate their affected eqs)
E=set(F0)
for knob in [7068,4432]:
    ae=sorted(affected_eqs(knob))
    touch=sorted(set(desc_of[knob]))
    r0={i:eval(H.eqcode[i],ns) for i in ae}
    setmany({knob:base[knob]+1}, touch)
    r1={i:eval(H.eqcode[i],ns) for i in ae}
    setmany({knob:base[knob]}, touch)
    broke=[i for i in ae if r1[i]!=r0[i] and r0[i]==0]
    E.update(broke)
E=sorted(E)
K=sorted(set(v for i in E for v in H.eqvars[i] if v in H.freeinp))
print(f"E={len(E)} equations, K={len(K)} free vars")
# affinity test over Z (only eval E), 2nd difference per var
touchall={v:sorted(set(desc_of[v])) for v in K}
def residE():
    return {i:eval(H.eqcode[i],ns) for i in E}
r0=residE()
nonlin=set()
for v in K:
    o=base[v]
    setmany({v:o+1},touchall[v]); r1=residE()
    setmany({v:o+2},touchall[v]); r2=residE()
    setmany({v:o},touchall[v])
    for i in E:
        if (r2[i]-r1[i])!=(r1[i]-r0[i]): nonlin.add(v); break
print(f"nonlinear vars in K: {len(nonlin)}: {sorted(nonlin)}")
lin=[v for v in K if v not in nonlin]
print(f"linear vars: {len(lin)}")
# joint affinity check among linear vars at random point
random.seed(5)
pt={v:base[v]+random.randint(-3,3) for v in lin}
tset=sorted(set(k for v in lin for k in touchall[v]))
for v in lin: H.val[v]=pt[v]
for k in tset: H.val[H.order[k]]=eval(H.gcode[k],ns)
ract=residE()
for v in lin: H.val[v]=base[v]
for k in tset: H.val[H.order[k]]=eval(H.gcode[k],ns)
# build slopes
slope={}
for v in lin:
    o=base[v]; setmany({v:o+1},touchall[v]); r1=residE(); setmany({v:o},touchall[v])
    slope[v]={i:r1[i]-r0[i] for i in E}
jointok=True
for i in E:
    pred=r0[i]+sum(slope[v][i]*(pt[v]-base[v]) for v in lin)
    if pred!=ract[i]: jointok=False;break
print(f"joint affine over linear vars exact? {jointok}")
json.dump({'E':E,'lin':lin,'nonlin':sorted(nonlin)},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/closure2.json','w'))
