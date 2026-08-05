import heal_harness as H
import json, pickle, time
from collections import defaultdict
p=H.p
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
ATOMS=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']])
d=H.loadd('best/new_instance_partial_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
sol=pickle.load(open(SCR+'/cert22.pkl','rb'))
breakable=sol['breakable']; nz=sol['nz']
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
def incr(w,nv):
    H.val[w]=nv
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],{'v':H.val,'__builtins__':{}})
def av(items):
    s=0
    for m,c in items:
        tt=c%p
        for v in m: tt=tt*H.val[v]%p
        s=(s+tt)%p
    return s
def free_anc_atom(i):
    s=set()
    for m,c in ATOMS[i]:
        for v in m:
            if v in H.freeinp: s.add(v)
            else: s|=H.anc.get(v,set())
    return s
# classify affine vs nonlinear (test 2nd diff on each free ancestor; cap for speed)
t0=time.time()
affine=set(); nonlin=set()
for a in breakable:
    fa=sorted(free_anc_atom(a))
    isl=True
    for f in fa:
        b0=av(ATOMS[a]); incr(f,base[f]+1); b1=av(ATOMS[a]); incr(f,base[f]+2); b2=av(ATOMS[a]); incr(f,base[f])
        if (b1-b0)%p!=(b2-b1)%p: isl=False; break
    (affine if isl else nonlin).add(a)
print(f"classified {len(breakable)} atoms in {time.time()-t0:.1f}s: affine={len(affine)}, nonlinear={len(nonlin)}")
pickle.dump({'affine':sorted(affine),'nonlin':sorted(nonlin)},open(SCR+'/classify22.pkl','wb'))

# Build Jacobian with ONLY affine atoms as constraints; target G2=atom20864 (and G1=20862).
# Reuse forward-mode from pipeline
gatepoly=pickle.load(open(SCR+'/gatepoly.pkl','rb'))
gate_items={t:list(pl.items()) for t,pl in gatepoly.items()}
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
for v in H.freeinp: H.val[v]=base[v]
H.forward(); val=H.val
localJac={}
for t,items in gate_items.items():
    localJac[t]={v:ppi(items,v,val) for v in set(x for m,c in items for x in m)}
# constraint atoms = affine breakable atoms; targets = G1,G2 (they are affine? check)
print("G1(20862) affine?",20862 in affine," G2(20864) affine?",20864 in affine)
useatoms=sorted(affine)  # includes G1,G2 if affine
atomJac={}; var2atom=defaultdict(list)
for a in useatoms:
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
r={a:av(ATOMS[a]) for a in useatoms}
# GE: wiring affine atoms (rhs0) first, then G1,G2 targets
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
            if rhs%p!=0: return False,(ri,combo,rhs),piv,ech
            continue
        pcol=min(row); inv=pow(row[pcol],p-2,p)
        piv[pcol]=len(ech); ech.append([{c:cf*inv%p for c,cf in row.items()},rhs*inv%p,{k:mv*inv%p for k,mv in combo.items()}])
    return True,None,piv,ech
targets=[a for a in [20862,20864] if a in useatoms]
wiring=[a for a in useatoms if Jrows.get(a) and r[a]==0 and a not in targets]
rows=[(Jrows[a],0) for a in wiring]+[(Jrows.get(a,{}),(-r[a])%p) for a in targets]
meta=wiring+targets
ok,inc,piv,ech=rref(rows)
print(f"AFFINE-ONLY solve: consistent={ok}, rank={len(piv)}, #affine-wiring={len(wiring)}, targets={targets}")
if not ok:
    ri,combo,rhs=inc
    print(f"  INCONSISTENT at atom {meta[ri]} -> G2/G1 AFFINELY PINNED (global infeasibility signal)")
    print(f"  cert size={sum(1 for k,mv in combo.items() if mv%p)}")
else:
    print("  CONSISTENT with affine atoms only -> nonlinearity was essential -> escapable!")
