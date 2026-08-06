import heal_harness as H
import json, pickle
from collections import defaultdict
p=H.p
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
ATOMS=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']])
gatepoly=pickle.load(open(SCR+'/gatepoly.pkl','rb'))
gate_items={t:list(pl.items()) for t,pl in gatepoly.items()}
cl=pickle.load(open(SCR+'/classify22.pkl','rb')); affine=set(cl['affine'])
sol=pickle.load(open(SCR+'/cert22.pkl','rb')); breakable=set(sol['breakable'])
d=H.loadd('best/new_instance_partial_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); val=H.val
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
def boolean_var(i):
    vs=set(x for m,c in ATOMS[i] for x in m)
    if len(vs)!=1: return None
    v=next(iter(vs)); c1=c2=0
    for m,c in ATOMS[i]:
        if len(m)==2: c2+=c
        elif len(m)==1: c1+=c
        elif len(m)==0 and c!=0: return None
    if c2!=0 and c1!=0 and (c1+c2)%p==0: return v
    return None
BOOLVARS=set(bv for i in range(len(ATOMS)) if (bv:=boolean_var(i)) is not None)
BOOLFREE=BOOLVARS&H.freeinp
def is_bool_atom(i):
    for m,c in ATOMS[i]:
        if len(m)==2 and m[0]!=m[1]: return False
    return any(len(m)==2 for m,c in ATOMS[i])
boolatoms=set(i for i in breakable if is_bool_atom(i))
def ppi(items,v):
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
def av(i):
    s=0
    for m,c in ATOMS[i]:
        tt=c%p
        for v in m: tt=tt*val[v]%p
        s=(s+tt)%p
    return s
localJac={}
for t,items in gate_items.items():
    localJac[t]={v:ppi(items,v) for v in set(x for m,c in items for x in m)}
KNOBS=[f for f in H.freeinp if f not in BOOLFREE]
def buildJ(useatoms):
    atomJac={}; var2atom=defaultdict(list)
    for a in useatoms:
        jac={}
        for v in set(x for m,c in ATOMS[a] for x in m):
            pv=ppi(ATOMS[a],v)
            if pv: jac[v]=pv
        atomJac[a]=jac
        for v in jac: var2atom[v].append(a)
    Jrows=defaultdict(dict)
    for f in KNOBS:
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
def rref(rows):
    piv={}; ech=[]
    for ri,(row0,rhs0) in enumerate(rows):
        row=dict(row0); rhs=rhs0%p; combo={ri:1}
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
            for k,mv in e[2].items(): combo[k]=(combo.get(k,0)-factor*mv)%p
        if not row:
            if rhs%p!=0: return False,(ri,combo)
            continue
        pcol=min(row); inv=pow(row[pcol],p-2,p)
        piv[pcol]=len(ech); ech.append([{c:cf*inv%p for c,cf in row.items()},rhs*inv%p,{k:mv*inv%p for k,mv in combo.items()}])
    return True,None
gap=[20862,20864]
loadatoms={18081,18084,29377,35321}
for label,excl in [("exclude booleans + atom18081", boolatoms|{18081}),
                   ("exclude booleans + 4 load atoms", boolatoms|loadatoms),
                   ("exclude booleans + load atoms + ripple(7450,7452)", boolatoms|loadatoms|{7450,7452})]:
    useatoms=sorted(breakable-excl)
    Jrows=buildJ(set(useatoms)|set(gap))
    r={a:av(a) for a in set(useatoms)|set(gap)}
    wiring=[a for a in useatoms if Jrows.get(a) and r[a]==0]
    rows=[(Jrows[a],0) for a in wiring]+[(Jrows.get(a,{}),(-r[a])%p) for a in gap]
    meta=wiring+gap
    ok,inc=rref(rows)
    msg=f"{label}: consistent={ok}"
    if not ok:
        ri,combo=inc; cert=[meta[k] for k,mv in combo.items() if mv%p]
        msg+=f"  cert size={len(cert)}, affine-only={all(a in affine for a in cert)}"
    print(msg)
