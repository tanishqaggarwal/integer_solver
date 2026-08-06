"""Exact forward-mode AD Jacobian of check atoms wrt free inputs at current H.val state."""
import heal_harness as H
import pickle, time
from collections import defaultdict
p=H.p
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'

gatepoly=pickle.load(open(SCR+'/gatepoly.pkl','rb'))   # target -> {monomial:coef}
ai=pickle.load(open(SCR+'/atoms.pkl','rb'))
atoms=ai['atoms']; breakable=ai['breakable']; nz=ai['nz']

# gate poly as list items
gate_items={t:list(pl.items()) for t,pl in gatepoly.items()}

def poly_partial_items(items, v, val):
    s=0
    for m,c in items:
        cnt=m.count(v)
        if cnt==0: continue
        term=(c*cnt)%p
        seen=False
        for u in m:
            if u==v and not seen: seen=True; continue
            term=term*val[u]%p
        s=(s+term)%p
    return s

def build_jac(val):
    t0=time.time()
    # local gate jacobian
    localJac={}
    for t,items in gate_items.items():
        vs=set()
        for m,c in items: vs.update(m)
        localJac[t]={v:poly_partial_items(items,v,val) for v in vs}
    # atom jacobian (check atoms only)
    atomJac={}
    var2atom=defaultdict(list)
    for a in breakable:
        items=atoms[a]
        vs=set()
        for m,c in items: vs.update(m)
        jac={v:poly_partial_items(items,v,val) for v in vs}
        # drop zero partials
        jac={v:c for v,c in jac.items() if c}
        atomJac[a]=jac
        for v in jac: var2atom[v].append(a)
    # desc_of per free input, topo order
    desc_of=defaultdict(list)
    for k,t in enumerate(H.order):
        for w in H.anc[t]:
            desc_of[w].append(k)
    # forward-mode per free input
    Jrows=defaultdict(dict)  # atom -> {free: coef}
    for f in H.freeinp:
        dof=desc_of.get(f)
        if not dof:
            # f may still directly appear in atoms as a var
            touched_atoms=set(var2atom.get(f,()))
            d={f:1}
        else:
            d={f:1}
            for k in dof:
                t=H.order[k]
                lj=localJac[t]
                acc=0
                for v,cf in lj.items():
                    dv=d.get(v)
                    if dv: acc=(acc+cf*dv)%p
                if acc: d[t]=acc
            touched_atoms=set()
            for v in d: touched_atoms.update(var2atom.get(v,()))
        for a in touched_atoms:
            jac=atomJac[a]
            acc=0
            for v,cf in jac.items():
                dv=d.get(v)
                if dv: acc=(acc+cf*dv)%p
            if acc: Jrows[a][f]=acc
    # residuals
    def atomval(items):
        s=0
        for m,c in items:
            tt=c%p
            for v in m: tt=tt*val[v]%p
            s=(s+tt)%p
        return s
    r={a:atomval(atoms[a]) for a in breakable}
    print(f"  jac built in {time.time()-t0:.1f}s; rows(atoms w/ nonzero jac)={sum(1 for a in Jrows if Jrows[a])}; nnz={sum(len(row) for row in Jrows.values())}")
    return Jrows, r

if __name__=='__main__':
    d=H.loadd('best/new_instance_partial_39013.json')
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward()
    Jrows,r=build_jac(H.val)
    nzr=[a for a in breakable if r[a]!=0]
    print("nonzero-residual atoms:", nzr)
    for a in nzr:
        print(f"  atom {a}: |jac row|={len(Jrows.get(a,{}))}, residual={r[a]}")
    pickle.dump({'Jrows':{a:dict(row) for a,row in Jrows.items()},'r':r,'breakable':breakable,'nz':nz},
                open(SCR+'/jac.pkl','wb'))
    print("saved jac.pkl")
