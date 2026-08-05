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
cl=pickle.load(open(SCR+'/classify22.pkl','rb')); affine=set(cl['affine']); nonlin=set(cl['nonlin'])
sol=pickle.load(open(SCR+'/cert22.pkl','rb')); breakable=set(sol['breakable'])
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
# identify boolean atoms: an atom is "boolean-type" if EVERY quadratic monomial (u,u) in it
# and it reduces to booleanness. Heuristic: atom repr contains only patterns x*(x-1), x^2-x, 2x(1-x)
# Structural: all monomials have len<=2, and for every degree-2 monomial (a,b) we have a==b (pure squares),
# i.e., no cross-products (bilinear). Bilinear (a,b),a!=b => genuine product (pin/verifier), keep.
def is_boolean_atom(i):
    for m,c in ATOMS[i]:
        if len(m)==2 and m[0]!=m[1]:
            return False  # has a cross product -> not pure-boolean
    # must have at least one square term
    return any(len(m)==2 for m,c in ATOMS[i])
boolatoms=set(i for i in breakable if is_boolean_atom(i))
print("boolean-type breakable atoms:", len(boolatoms))
# non-boolean nonlinear = nonlin minus boolatoms
nonbool_nonlin=nonlin-boolatoms
print("non-boolean nonlinear breakable atoms:", len(nonbool_nonlin))
d=H.loadd('best/new_instance_partial_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); val=H.val
localJac={}
for t,items in gate_items.items():
    localJac[t]={v:ppi(items,v,val) for v in set(x for m,c in items for x in m)}
def jacrow_for(atomset):
    atomJac={}; var2atom=defaultdict(list)
    for a in atomset:
        jac={}
        for v in set(x for m,c in ATOMS[a] for x in m):
            pv=ppi(ATOMS[a],v,val)
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
    return Jrows
def consistent(rows, track=False):
    piv={}; ech=[]
    for ri,(row0,rhs0) in enumerate(rows):
        row=dict(row0); rhs=rhs0%p; combo={ri:1} if track else None
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
            if track:
                for k,mv in e[2].items(): combo[k]=(combo.get(k,0)-factor*mv)%p
        if not row:
            if rhs%p!=0: return False,len(piv),(ri,combo)
            continue
        pcol=min(row); inv=pow(row[pcol],p-2,p)
        piv[pcol]=len(ech); ech.append([{c:cf*inv%p for c,cf in row.items()},rhs*inv%p, {k:mv*inv%p for k,mv in combo.items()} if track else None])
    return True,len(piv),None
gap=[20862,20864]
# constraint set = affine + non-boolean nonlinear (exclude booleans)
keepset=sorted((affine|nonbool_nonlin))
Jrows=jacrow_for(set(keepset)|set(gap))
r={a:av(ATOMS[a],val) for a in set(keepset)|set(gap)}
wiring=[a for a in keepset if Jrows.get(a) and r[a]==0]
rows=[(Jrows[a],0) for a in wiring]+[(Jrows.get(a,{}),(-r[a])%p) for a in gap]
meta=wiring+gap
ok,rank,inc=consistent([(row,rhs) for row,rhs in rows],track=True)
print(f"\nEXCLUDE booleans: consistent={ok}, rank={rank}, #constraints={len(wiring)}")
if not ok:
    ri,combo=inc
    cert=[(meta[k],mv%p) for k,mv in combo.items() if mv%p]
    print("  still inconsistent; cert atoms:",[a for a,_ in cert][:30])
    # classify cert atoms
    print("  cert affine:",sum(1 for a,_ in cert if a in affine),"nonlin:",sum(1 for a,_ in cert if a in nonlin))
    pickle.dump({'cert':cert},open(SCR+'/cert_nobool.pkl','wb'))
