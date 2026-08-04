import heal_harness as H
from propagate import load_atoms, atom_vars
from collections import defaultdict
p=H.p
A=load_atoms()
mem=defaultdict(int)
for poly in A:
    for v in atom_vars(poly): mem[v]+=1
def evalp(poly,V):
    s=0
    for m,c in poly.items():
        t=c
        for v in m: t*=V[v]
        s+=t
    return s
v013=H.loadd('best/new_instance_partial_39013.json'); vA=H.loadd('best_agentA_39022.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
V=H.val[:]
baseF=set(H.fails())
ns={'v':H.val,'__builtins__':{}}
# collect all product-slack partner candidates from a44342, a45677
# a candidate knob = a var v s.t. it appears in a product term x_v*x_w where x_w != 0 (so linear effect)
cands=set()
for ai in [44342,45677]:
    for m,c in A[ai].items():
        if len(m)==2:
            a,b=m
            if V[a]!=0 and V[b]==0: cands.add(b)
            if V[b]!=0 and V[a]==0: cands.add(a)
        if len(m)==1:
            if V[m[0]]==0: cands.add(m[0])
print(f"linear-effect knob candidates in a44342/a45677: {sorted(cands)}")
# for each, compute TRUE fanout: perturb by 1, forward, count broken satisfied eqs
def fanout(v):
    old=V_full = None
    H.val[v]+=1
    H.forward()
    newF=set(H.fails())
    broke=len(newF-baseF)
    fixedcount=len(baseF-newF)
    # restore
    for k in range(len(H.val)): H.val[k]=V[k]
    H.forward()
    return broke, sorted(newF-baseF)[:8]
print("var    free  atoms  fanout(new-broken)  which")
rows=[]
for v in sorted(cands):
    b,wh=fanout(v)
    rows.append((b,v,wh))
rows.sort()
for b,v,wh in rows:
    print(f"x_{v:<6} {str(v in H.freeinp):<5} {mem[v]:<5} broke={b:<4} {wh}")
