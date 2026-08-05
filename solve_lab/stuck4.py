import heal_harness as H
import json, pickle
from collections import defaultdict
p=H.p
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
ATOMS=[]; reprs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr',''))
gatepoly=pickle.load(open(SCR+'/gatepoly.pkl','rb'))
gate_items={t:list(pl.items()) for t,pl in gatepoly.items()}
cl=pickle.load(open(SCR+'/classify22.pkl','rb')); affine=set(cl['affine'])
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
def ppi(items,v,val):
    s=0
    for m,c in items:
        cnt=m.count(v)
        if cnt==0: continue
        term=(c*cnt)%p; seen=False
        for u in m:
            if u==v and not seen: seen=True; continue
            term=term*val[u]%p
        s=(s+term)%p
    return s
def av(items,val):
    s=0
    for m,c in items:
        tt=c%p
        for v in m: tt=tt*val[v]%p
        s=(s+tt)%p
    return s
def solve_particular(rows):
    piv={}; ech=[]
    for (row0,rhs0) in rows:
        row=dict(row0); rhs=rhs0%p
        while True:
            pc=None
            for c in row:
                if c in piv: pc=c;break
            if pc is None: break
            e=ech[piv[pc]]; factor=row[pc]%p
            for c,cf in e[0].items():
                nv=(row.get(c,0)-factor*cf)%p
                if nv: row[c]=nv
                elif c in row: del row[c]
            rhs=(rhs-factor*e[1])%p
        if not row:
            if rhs%p!=0: return None
            continue
        pcol=min(row); inv=pow(row[pcol],p-2,p)
        piv[pcol]=len(ech); ech.append([{c:cf*inv%p for c,cf in row.items()},rhs*inv%p])
    x={}
    for _ in range(len(ech)+3):
        ch=False
        for pcol,idx in piv.items():
            row,rhs=ech[idx]; s=rhs
            for c,cf in row.items():
                if c==pcol: continue
                s=(s-cf*x.get(c,0))%p
            if x.get(pcol)!=s%p: x[pcol]=s%p; ch=True
        if not ch: break
    return x

# redo one affine step to reach the stuck state
d=H.loadd('best/new_instance_partial_39022.json')
base=[0]*H.NVARS
for v in H.freeinp: base[v]=d.get(v,0)
H.val[:]=base; H.forward()
useatoms=sorted(affine)
localJac={}
for t,items in gate_items.items():
    localJac[t]={v:ppi(items,v,H.val) for v in set(x for m,c in items for x in m)}
atomJac={}; var2atom=defaultdict(list)
for a in useatoms:
    jac={}
    for v in set(x for m,c in ATOMS[a] for x in m):
        pv=ppi(ATOMS[a],v,H.val)
        if pv: jac[v]=pv
    atomJac[a]=jac
    for v in jac: var2atom[v].append(a)
Jrows=defaultdict(dict)
for f in H.freeinp:
    dof=desc_of.get(f); dd_={f:1}
    if dof:
        for k in dof:
            t=H.order[k]; acc=0
            for v,cf in localJac[t].items():
                dv=dd_.get(v)
                if dv: acc=(acc+cf*dv)%p
            if acc: dd_[t]=acc
    touched=set()
    for v in dd_: touched.update(var2atom.get(v,()))
    for a in touched:
        acc=0
        for v,cf in atomJac[a].items():
            dv=dd_.get(v)
            if dv: acc=(acc+cf*dv)%p
        if acc: Jrows[a][f]=acc
r={a:av(ATOMS[a],H.val) for a in useatoms}
targets=[a for a in useatoms if r[a]!=0]
wiring=[a for a in useatoms if Jrows.get(a) and r[a]==0]
rows=[(Jrows[a],0) for a in wiring]+[(Jrows.get(a,{}),(-r[a])%p) for a in targets]
x=solve_particular(rows)
for f,dv in x.items():
    if f in H.freeinp: H.val[f]=(H.val[f]+dv)%p
H.forward()
nzall=[a for a in range(len(ATOMS)) if av(ATOMS[a],H.val)!=0]
print("stuck nonzero atoms:",nzall)
for a in nzall:
    vs=sorted(set(x for m,c in ATOMS[a] for x in m))
    print(f"\n--- atom {a} ({'affine' if a in affine else 'NONLINEAR'}): {reprs[a][:120]}")
    for v in vs[:14]:
        print(f"    x_{v}={H.val[v]%p}  free={v in H.freeinp} anc={sorted(H.anc.get(v,set()))[:5]}")
