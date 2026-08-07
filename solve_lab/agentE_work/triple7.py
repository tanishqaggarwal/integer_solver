"""Simultaneous exact solve of the triple AND its collateral, over the AFFINE knobs only."""
import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(5000000)
import engine as E, fast, sparse, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
s={int(k):int(v) for k,v in json.load(open('triple_state_seed.json')).items()}
v0=E.forward(s); bad0=E.badatoms(v0)
print("start bad",sorted(bad0),flush=True)
def cone_free(a): return set(E.cone(a)[1])
cand=set()
for a in bad0: cand|=cone_free(a)
# plus handles and the knob cones
for u in (25848,17317,18682,28841,27928,9776,3271,31154,6635,36280,19097):
    c=set(); st=[u]; seen=set()
    while st:
        w=st.pop()
        if w in seen: continue
        seen.add(w); dv=E.definer[w]
        if dv is None: c.add(w); continue
        for z in E.avars[dv[0]]:
            if z!=w: st.append(z)
    cand|=c
cand-= {18956,1530,1603}
cand=sorted(cand)
print("candidate knobs:",len(cand),flush=True)
cols={}; affine=[]
t0=time.time()
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
        if b2.get(a,0)-bad0.get(a,0)!=2*col.get(a,0): ok=False; break
        if b7.get(a,0)-bad0.get(a,0)!=7*col.get(a,0): ok=False; break
    if ok:
        affine.append(f); cols[f]=col
print(f"affine knobs: {len(affine)} of {len(cand)}  ({time.time()-t0:.0f}s)",flush=True)
atoms=set(bad0)
for f in affine: atoms|=set(cols[f])
atoms=sorted(atoms)
print("atoms in the affine system:",len(atoms),atoms[:30],flush=True)
rowmap={a:{} for a in atoms}
for f in affine:
    for a,c in cols[f].items(): rowmap[a][f]=c
rows=[rowmap[a] for a in atoms]; rhs=[-bad0.get(a,0) for a in atoms]
sol,msg,_=sparse.solve_sparse(rows,rhs,names=atoms,verbose=True,maxcore=400)
print("AFFINE SYSTEM ->",msg,flush=True)
if sol is not None:
    ns=dict(s)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"EXACT: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)}",flush=True)
    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('triple7_%d.json'%(39033-len(ff)),'w'))
    json.dump({str(k):str(int(x)) for k,x in ns.items()}, open('triple7_seed.json','w'))
