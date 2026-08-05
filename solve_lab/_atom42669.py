import heal_harness as H, json, random
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
a=atoms[42669]
print("atom 42669 full repr:")
print(a['repr'])
print("poly monomials:",a['poly'])
print("eqs:",a['eqs'])
# verify x_642 == x_28599 * x_17325 as gate identity (over random assignments)
gates={}
for line in open('atoms/gates.jsonl'):
    g=json.loads(line); gates[g['t']]=(g['rhs'],tuple(g['vids']))
print("\nx_642 gate:",gates.get(642))
print("x_28730 gate:",gates.get(28730))
# Test: are ALL 11 failing equation residuals exactly linear in free inputs (branch A)?
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
base={v:d.get(v,0) for v in H.freeinp}
selectors={4287,2081,9062,20434}
perturbable=[v for v in sorted(H.freeinp) if v not in selectors]
def setall(assign):
    for v in H.freeinp: H.val[v]=assign[v]
    H.forward()
def residF():
    ns={'v':H.val,'__builtins__':{}}
    return [eval(H.eqcode[i],ns)%p for i in F]
random.seed(5); nonlin=0
for trial in range(30):
    da={v:random.randint(-9,9) for v in random.sample(perturbable,60)}
    db={v:random.randint(-9,9) for v in random.sample(perturbable,60)}
    A=dict(base);
    for v,x in da.items(): A[v]=base[v]+x
    setall(A); ra=residF()
    B=dict(base)
    for v,x in db.items(): B[v]=base[v]+x
    setall(B); rb=residF()
    AB=dict(base)
    for v,x in da.items(): AB[v]=base[v]+x
    for v,x in db.items(): AB[v]=AB[v]+x
    setall(AB); rab=residF()
    setall(base); r0=residF()
    for k in range(len(F)):
        if (rab[k]-ra[k]-rb[k]+r0[k])%p!=0: nonlin+=1; break
print(f"\n11-fail residuals: nonlinear second-difference in {nonlin}/30 branch-A trials (0 => ALL 11 fails EXACTLY linear in free inputs)")
