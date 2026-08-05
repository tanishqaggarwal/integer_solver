import heal_harness as H
import pickle, json
p=H.p
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
data=pickle.load(open(SCR+'/jac.pkl','rb'))
Jrows=data['Jrows']; r=data['r']; breakable=data['breakable']
cd=pickle.load(open(SCR+'/cert.pkl','rb'))
combo=cd['combo']; rowmeta=cd['rowmeta']; cert=cd['cert']

# Reconstruct the exact row list used in gf_solve (must match ordering!)
resid_atoms=[a for a in breakable if r.get(a,0)!=0]
wiring=[a for a in breakable if a in Jrows and Jrows[a] and r.get(a,0)==0]
rows=[]
for a in wiring: rows.append((Jrows[a],0))
for a in resid_atoms: rows.append((Jrows.get(a,{}),(-r[a])%p))

# Verify: sum combo[k]*rows[k].LHS == 0 ; sum combo[k]*rows[k].RHS == cert rhs
lhs={}
rhs=0
for k,mv in combo.items():
    row,b=rows[k]
    for c,cf in row.items():
        lhs[c]=(lhs.get(c,0)+mv*cf)%p
    rhs=(rhs+mv*b)%p
lhs={c:v for c,v in lhs.items() if v%p}
print("certificate LHS nonzero cols:", len(lhs))
print("certificate RHS (should be nonzero):", rhs%p)
print("combo[5491] (atom18081 coef):", combo.get(5491))

# which atoms in certificate, and their degree in free inputs
cert_atoms=[(rowmeta[k], mv%p) for k,mv in combo.items() if mv%p]
print("num atoms in cert:", len(cert_atoms))
# load atom reprs + polys
reprs={}; apoly={}
with open('atoms/poly_atoms.jsonl') as f:
    for i,line in enumerate(f):
        dd=json.loads(line)
        reprs[i]=dd.get('repr','')
        apoly[i]=[(tuple(m),c) for m,c in dd['poly']]

# degree of each cert atom in FREE inputs: check if atom residual is linear in free.
# We test: is the atom's dependence on free inputs affine? Use finite diff order-2 on a few free.
d=H.loadd('best/new_instance_partial_39013.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
from collections import defaultdict
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
def incr(w,nv):
    H.val[w]=nv
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],{'v':H.val,'__builtins__':{}})
def restore():
    for v in H.freeinp: H.val[v]=base[v]
    H.forward()
def atomval(i):
    s=0
    for m,c in apoly[i]:
        tt=c%p
        for v in m: tt=tt*H.val[v]%p
        s=(s+tt)%p
    return s
def free_anc_of_atom(i):
    s=set()
    for m,c in apoly[i]:
        for v in m:
            if v in H.freeinp: s.add(v)
            else: s|=H.anc.get(v,set())
    return s

# classify each cert atom as affine or nonlinear in free inputs (test all its free ancestors? sample)
import random
random.seed(0)
lin_atoms=0; nonlin_atoms=0; nonlin_list=[]
for a,mv in cert_atoms:
    fa=sorted(free_anc_of_atom(a))
    islin=True
    for f in fa:
        b0=atomval(a)
        incr(f,base[f]+1); b1=atomval(a)
        incr(f,base[f]+2); b2=atomval(a)
        incr(f,base[f])  # restore this f (others unchanged since we only moved f)
        if (b1-b0)%p != (b2-b1)%p:
            islin=False; break
    if islin: lin_atoms+=1
    else: nonlin_atoms+=1; nonlin_list.append(a)
print(f"cert atoms affine-in-free: {lin_atoms}, nonlinear: {nonlin_atoms}")
print("nonlinear cert atoms:", nonlin_list[:20])
for a in nonlin_list[:8]:
    print(f"   atom {a}: {reprs[a]}")
print()
print("=== sample cert atom reprs ===")
for a,mv in cert_atoms[:25]:
    print(f" atom {a} (mult {mv%p if (mv%p)<10**6 else hex(mv%p)[:12]+'...'}): {reprs[a]}")
