import heal_harness as H, re, random, json
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
def setmany(assign):
    for v,x in assign.items(): H.val[v]=x
    ks=set()
    for v in assign: ks.update(desc_of[v])
    for k in sorted(ks): H.val[H.order[k]]=eval(H.gcode[k],ns)
def resid(idxs):
    return {i:eval(H.eqcode[i],ns) for i in idxs}
F0=set(H.fails())
# closure: find core free vars + affected equations
K=set([7068,4432,17325,9413])
allE=set(F0)
for _ in range(6):
    newE=set()
    for v in list(K):
        old=base[v]; setmany({v:old+1}); f1=set(H.fails()); setmany({v:old})
        newE|= (f1 ^ F0)
    allE |= newE
    # free vars in allE
    newK=set()
    for i in allE:
        for v in H.eqvars[i]:
            if v in H.freeinp: newK.add(v)
    if newK<=K and newE<=allE: 
        K|=newK; break
    K|=newK
print(f"closure: {len(K)} core free vars, {len(allE)} affected equations")
K=sorted(K); allE=sorted(allE)
# test joint affinity over Z: random integer perturbation
r0=resid(allE)
random.seed(11)
NT=3
allaff=True; nonlin_vars=set()
# first per-variable second difference to find nonlinear vars
for v in K:
    o=base[v]
    setmany({v:o+1}); r1=resid(allE)
    setmany({v:o+2}); r2=resid(allE)
    setmany({v:o})
    for i in allE:
        if (r2[i]-r1[i])!=(r1[i]-r0[i]): nonlin_vars.add(v); break
print(f"nonlinear core vars: {len(nonlin_vars)}: {sorted(nonlin_vars)}")
print(f"linear core vars: {len(K)-len(nonlin_vars)}")
json.dump({'K':K,'E':allE,'nonlin':sorted(nonlin_vars)}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/closure.json','w'))
