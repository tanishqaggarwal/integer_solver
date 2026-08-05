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
d=H.loadd('best/new_instance_partial_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); val=H.val
localJac={}
for t,items in gate_items.items():
    localJac[t]={v:ppi(items,v,val) for v in set(x for m,c in items for x in m)}
# jac rows for a set of atoms
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
def consistent(rows):
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
            if rhs%p!=0: return False,len(piv)
            continue
        pcol=min(row); inv=pow(row[pcol],p-2,p)
        piv[pcol]=len(ech); ech.append([{c:cf*inv%p for c,cf in row.items()},rhs*inv%p])
    return True,len(piv)

loadatoms=[18081,18084,29377,35321]
gap=[20862,20864]
allset=sorted(affine|set(loadatoms))
Jrows=jacrow_for(allset)
r={a:av(ATOMS[a],val) for a in allset}
# Test A: affine wiring + gap targets (no load constraints)
def run(extra_keep):
    wiring=[a for a in affine if Jrows.get(a) and r[a]==0]
    rows=[(Jrows[a],0) for a in wiring]
    for a in extra_keep:  # keep these nonlinear atoms at 0 (linearized)
        rows.append((Jrows.get(a,{}),0))
    for a in gap:  # drive gap to 0
        rows.append((Jrows.get(a,{}),(-r[a])%p))
    return consistent(rows)
print("A: affine + gap targets (no load keep):", run([]))
print("B: affine + keep x_11150(18081) + gap:", run([18081]))
print("C: affine + keep all 4 load atoms + gap:", run(loadatoms))
# Also: can we zero BOTH gap AND loads if loads were nonzero? test keeping loads while moving gap
# Reverse: from a state where loads are the obstruction, does zeroing loads force gap?
