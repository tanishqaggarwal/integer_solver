import heal_harness as H, random
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
selectors={4287,2081,9062,20434}
perturbable=[v for v in sorted(H.freeinp) if v not in selectors]
NEQ=len(H.eqcode)
def setall(assign):
    for v in H.freeinp: H.val[v]=assign[v]
    H.forward()
def residall():
    ns={'v':H.val,'__builtins__':{}}
    return [eval(c,ns)%p for c in H.eqcode]
nonlin=set()
random.seed(0)
NT=12
for trial in range(NT):
    S=random.sample(perturbable, 200)
    da={v:random.randint(1,10**6) for v in S}
    db={v:random.randint(1,10**6) for v in S}
    A=dict(base)
    for v,x in da.items(): A[v]=base[v]+x
    setall(A); ra=residall()
    B=dict(base)
    for v,x in db.items(): B[v]=base[v]+x
    setall(B); rb=residall()
    AB=dict(base)
    for v,x in da.items(): AB[v]=base[v]+x
    for v,x in db.items(): AB[v]=AB[v]+x
    setall(AB); rab=residall()
    setall(base); r0=residall()
    for i in range(NEQ):
        if (rab[i]-ra[i]-rb[i]+r0[i])%p!=0: nonlin.add(i)
    print(f"trial {trial}: cumulative nonlinear eqs = {len(nonlin)}",flush=True)
setall(base)
print(f"\nTOTAL nonlinear (branch-A) equations detected: {len(nonlin)} / {NEQ}")
F=set([2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125])
print("any of the 11 fails nonlinear?", sorted(nonlin&F))
import json
json.dump(sorted(nonlin), open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/nonlin_eqs.json','w'))
print("saved nonlin_eqs.json")
