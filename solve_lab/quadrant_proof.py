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
# Identify boolean vars: atom is exactly a single-var boolean v(v-1)/v^2-v/2v(1-v)
def boolean_var(i):
    vs=set(x for m,c in ATOMS[i] for x in m)
    if len(vs)!=1: return None
    v=next(iter(vs))
    # degrees present
    c1=c2=0
    for m,c in ATOMS[i]:
        if len(m)==2: c2+=c
        elif len(m)==1: c1+=c
        # constant term must be 0
        elif len(m)==0 and c!=0: return None
    if c2!=0 and c1!=0 and (c1+c2)%p==0:  # v^2*c2+v*c1 with c1=-c2 => c2(v^2-v)
        return v
    return None
BOOLVARS=set()
for i in range(len(ATOMS)):
    bv=boolean_var(i)
    if bv is not None: BOOLVARS.add(bv)
BOOLFREE=BOOLVARS & H.freeinp
print("boolean vars:",len(BOOLVARS)," boolean FREE inputs (control bits):",len(BOOLFREE))
# boolean atoms (any atom whose only nonlinearity is boolean squares)
def is_bool_atom(i):
    for m,c in ATOMS[i]:
        if len(m)==2 and m[0]!=m[1]: return False
    return any(len(m)==2 for m,c in ATOMS[i])
boolatoms=set(i for i in breakable if is_bool_atom(i))

d=H.loadd('best/new_instance_partial_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); val=H.val
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
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
# FIX control bits: exclude BOOLFREE from the perturbation columns.
KNOBS=[f for f in H.freeinp if f not in BOOLFREE]
# constraint atoms: all breakable EXCEPT pure boolean atoms (auto-satisfied when bits fixed)
useatoms=sorted(breakable-boolatoms)
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
r={a:av(a) for a in useatoms}
gap=[20862,20864]
wiring=[a for a in useatoms if Jrows.get(a) and r[a]==0]
rows=[(Jrows[a],0) for a in wiring]+[(Jrows.get(a,{}),(-r[a])%p) for a in gap]
meta=wiring+gap
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
ok,inc=rref(rows)
print(f"\nFIX control bits + exclude boolean atoms: consistent={ok}")
if not ok:
    ri,combo=inc
    cert=[(meta[k],mv%p) for k,mv in combo.items() if mv%p]
    naff=sum(1 for a,_ in cert if a in affine); nnl=sum(1 for a,_ in cert if a not in affine)
    print(f"  INCONSISTENT. certificate: {len(cert)} atoms, affine={naff}, nonlinear-classified={nnl}")
    # rigorously verify EACH cert atom is affine in the KNOBS (non-bit free inputs)
    nonbitfree=set(KNOBS)
    def free_anc_atom(i):
        s=set()
        for m,c in ATOMS[i]:
            for v in m:
                if v in H.freeinp: s.add(v)
                else: s|=H.anc.get(v,set())
        return s
    def incr(w,nv):
        val[w]=nv
        for k in desc_of[w]: val[H.order[k]]=eval(H.gcode[k],{'v':val,'__builtins__':{}})
    base2={v:val[v] for v in H.freeinp}
    allaff=True
    for a,mv in cert:
        fa=[f for f in free_anc_atom(a) if f in nonbitfree]  # only non-bit knobs
        for f in fa:
            b0=av(a); incr(f,base2[f]+1); b1=av(a); incr(f,base2[f]+2); b2=av(a); incr(f,base2[f])
            if (b1-b0)%p!=(b2-b1)%p: allaff=False; print("   cert atom",a,"nonlinear in knob",f); break
    print(f"  ==> ALL cert atoms affine in non-bit knobs: {allaff}")
    if allaff:
        print("  ==> RIGID AFFINE CERTIFICATE: this control-bit quadrant is PROVABLY mod-p INFEASIBLE.")
    pickle.dump({'cert':cert,'meta':meta},open(SCR+'/cert_quadrant.pkl','wb'))
