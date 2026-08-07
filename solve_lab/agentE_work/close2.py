"""Close a20215 and a28647 by an exact affine solve over the free variables in the cones of
   x_24908 and of (x_6083, x_33708), with NON-BOOLEAN integer handles explicitly included."""
import sys, json, time, re
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, sparse, harness as H

SEED=sys.argv[1] if len(sys.argv)>1 else 'triple8_seed.json'
s={int(k):int(v) for k,v in json.load(open(SEED)).items()}
v0=E.forward(s); bad0=E.badatoms(v0)
print("start bad",sorted(bad0),"fails",len(E.eqfails(bad0)),flush=True)

def varcone(u):
    c=set(); st=[u]; seen=set()
    while st:
        w=st.pop()
        if w in seen: continue
        seen.add(w); dv=E.definer[w]
        if dv is None: c.add(w); continue
        for z in E.avars[dv[0]]:
            if z!=w: st.append(z)
    return c,seen

def is_bool(f):
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in ('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)'): return True
    return False

cand=set()
for u in (24908, 6083, 33708, 24530, 36433, 36990, 19239, 26386, 27475, 5647):
    c,_=varcone(u); cand|=c
for a in bad0: cand|=set(E.cone(a)[1])
cand-= {18956,1530,1603}
cand=sorted(cand)
nb=[f for f in cand if not is_bool(f)]
print(f"candidate knobs {len(cand)}  (non-boolean {len(nb)}, boolean {len(cand)-len(nb)})",flush=True)

t0=time.time(); cols={}; affine=[]; nonaff=0
for f in cand:
    o=v0[f]
    b1,_=fast.resid_delta(v0,bad0,{f:o+1})
    b2,_=fast.resid_delta(v0,bad0,{f:o+2})
    b7,_=fast.resid_delta(v0,bad0,{f:o+7})
    col={}
    for a in set(b1)|set(bad0):
        d=b1.get(a,0)-bad0.get(a,0)
        if d: col[a]=d
    ok=True
    for a in set(b2)|set(b7)|set(col)|set(bad0):
        if b2.get(a,0)-bad0.get(a,0)!=2*col.get(a,0) or b7.get(a,0)-bad0.get(a,0)!=7*col.get(a,0):
            ok=False; break
    if ok: affine.append(f); cols[f]=col
    else: nonaff+=1
print(f"affine knobs {len(affine)} of {len(cand)} ({time.time()-t0:.0f}s); "
      f"non-boolean among them: {sum(1 for f in affine if not is_bool(f))}",flush=True)
for a in (20215,28647):
    movers=[f for f in affine if a in cols[f]]
    print(f"  affine knobs moving a{a}: {len(movers)} -> {movers[:25]}",flush=True)

atoms=sorted(set(bad0)|set().union(*[set(cols[f]) for f in affine]) if affine else set(bad0))
print("atoms in the affine system:",len(atoms),flush=True)
rowmap={a:{} for a in atoms}
for f in affine:
    for a,c in cols[f].items(): rowmap[a][f]=c
rows=[rowmap[a] for a in atoms]; rhs=[-bad0.get(a,0) for a in atoms]
sol,msg,_=sparse.solve_sparse(rows,rhs,names=atoms,verbose=True,maxcore=600)
print("FULL AFFINE SYSTEM ->",msg,flush=True)
if sol is None:
    keep=[]
    for i,a in enumerate(atoms):
        idx=keep+[i]
        s2,m2,_=sparse.solve_sparse([rowmap[atoms[j]] for j in idx],[-bad0.get(atoms[j],0) for j in idx],
                                    verbose=False,maxcore=600)
        if s2 is not None: keep=idx
    drop=[atoms[i] for i in range(len(atoms)) if i not in keep]
    print("max solvable rows",len(keep),"of",len(atoms),"  BLOCKING:",drop,flush=True)
    for a in drop: print("   ATOM",a,H.atoms[a][:170],flush=True)
    sol,msg,_=sparse.solve_sparse([rowmap[atoms[j]] for j in keep],[-bad0.get(atoms[j],0) for j in keep],
                                  verbose=False,maxcore=600)
if sol is not None:
    ns=dict(s)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"EXACT: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)}",flush=True)
    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('close2_%d.json'%(39033-len(ff)),'w'))
    json.dump({str(k):str(int(x)) for k,x in ns.items()}, open('close2_seed.json','w'))
