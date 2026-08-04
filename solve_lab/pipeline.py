"""Full pipeline at a given state: build exact atom Jacobian, GE consistency, solve."""
import heal_harness as H
import pickle, time, json, sys
from collections import defaultdict
p=H.p
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
gatepoly=pickle.load(open(SCR+'/gatepoly.pkl','rb'))
gate_items={t:list(pl.items()) for t,pl in gatepoly.items()}
# load atoms
ATOMS=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']])
# breakable set (state-independent structure). Recompute once via random perturbation on a base.
def load_state(path):
    d=H.loadd(path)
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward()
def poly_partial_items(items,v,val):
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
def atomval(items,val):
    s=0
    for m,c in items:
        tt=c%p
        for v in m: tt=tt*val[v]%p
        s=(s+tt)%p
    return s
# desc_of (state-independent)
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)

def get_breakable():
    import random
    base=H.val[:]
    br=set()
    for seed in range(4):
        random.seed(500+seed)
        for v in H.freeinp: H.val[v]=random.randrange(p)
        H.forward()
        for i,a in enumerate(ATOMS):
            if atomval(a,H.val)!=0: br.add(i)
        H.val[:]=base; H.forward()
    return sorted(br)

def build_jac(breakable, val):
    localJac={}
    for t,items in gate_items.items():
        vs=set()
        for m,c in items: vs.update(m)
        localJac[t]={v:poly_partial_items(items,v,val) for v in vs}
    atomJac={}; var2atom=defaultdict(list)
    for a in breakable:
        jac={}
        for v in set(x for m,c in ATOMS[a] for x in m):
            pv=poly_partial_items(ATOMS[a],v,val)
            if pv: jac[v]=pv
        atomJac[a]=jac
        for v in jac: var2atom[v].append(a)
    Jrows=defaultdict(dict)
    for f in H.freeinp:
        dof=desc_of.get(f)
        d={f:1}
        if dof:
            for k in dof:
                t=H.order[k]; acc=0
                for v,cf in localJac[t].items():
                    dv=d.get(v)
                    if dv: acc=(acc+cf*dv)%p
                if acc: d[t]=acc
        touched=set()
        for v in d: touched.update(var2atom.get(v,()))
        for a in touched:
            acc=0
            for v,cf in atomJac[a].items():
                dv=d.get(v)
                if dv: acc=(acc+cf*dv)%p
            if acc: Jrows[a][f]=acc
    r={a:atomval(ATOMS[a],val) for a in breakable}
    return Jrows, r

def rref(rows, track=False):
    piv={}; ech=[]
    for ri,(row0,rhs0) in enumerate(rows):
        row=dict(row0); rhs=rhs0%p
        combo={ri:1} if track else None
        while True:
            pc=None
            for c in row:
                if c in piv: pc=c; break
            if pc is None: break
            e=ech[piv[pc]]; factor=row[pc]%p; er=e[0]
            for c,cf in er.items():
                nv=(row.get(c,0)-factor*cf)%p
                if nv: row[c]=nv
                elif c in row: del row[c]
            rhs=(rhs-factor*e[1])%p
            if track:
                for k,mv in e[2].items(): combo[k]=(combo.get(k,0)-factor*mv)%p
        if not row:
            if rhs%p!=0:
                return False,(ri,combo,rhs),piv,ech
            continue
        pcol=min(row); inv=pow(row[pcol],p-2,p)
        nrow={c:cf*inv%p for c,cf in row.items()}
        piv[pcol]=len(ech)
        ech.append([nrow,(rhs*inv)%p,{k:mv*inv%p for k,mv in combo.items()} if track else None])
    return True,None,piv,ech

if __name__=='__main__':
    path=sys.argv[1] if len(sys.argv)>1 else 'best/new_instance_partial_39022.json'
    load_state(path)
    print("state:",path,"fails:",len(H.fails()))
    breakable=get_breakable()
    print("breakable atoms:",len(breakable))
    nz=[a for a in breakable if atomval(ATOMS[a],H.val)!=0]
    print("nonzero atoms:",nz)
    t0=time.time()
    Jrows,r=build_jac(breakable,H.val)
    print(f"jac built {time.time()-t0:.1f}s; nnz={sum(len(v) for v in Jrows.values())}")
    # rows: wiring (rhs0) first, then nonzero-resid
    wiring=[a for a in breakable if Jrows.get(a) and r[a]==0]
    resid=[a for a in breakable if r[a]!=0]
    rows=[(Jrows[a],0) for a in wiring]+[(Jrows.get(a,{}),(-r[a])%p) for a in resid]
    meta=wiring+resid
    ok,inc,piv,ech=rref(rows,track=True)
    print(f"consistent={ok}, rank={len(piv)}")
    if not ok:
        ri,combo,rhs=inc
        print(f"INCONSISTENT at row {ri} = atom {meta[ri]}, rhs={rhs}")
        cert=[(meta[k],mv%p) for k,mv in combo.items() if mv%p]
        print("cert size:",len(cert))
        pickle.dump({'cert':cert,'meta':meta,'combo':combo,'inc':inc,'breakable':breakable,'nz':nz},open(SCR+'/cert22.pkl','wb'))
    else:
        pickle.dump({'piv':piv,'ech':ech,'meta':meta,'rows':rows,'breakable':breakable,'nz':nz,'r':r,'Jrows':{a:dict(v) for a,v in Jrows.items()}},open(SCR+'/solve22.pkl','wb'))
        print("saved solve22.pkl")
