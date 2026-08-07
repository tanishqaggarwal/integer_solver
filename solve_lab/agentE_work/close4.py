"""Joint exact solve of the 4-atom cluster {20215, 28647, 7389, 10187}: knobs are every free
   variable in their cones on which THOSE residuals are exactly affine (+1,+2,+7)."""
import sys, json, time, re
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, sparse, harness as H
CL=[20215,28647,7389,10187]
s={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}
v0=E.forward(s); bad0=E.badatoms(v0)
print("start bad",sorted(bad0),flush=True)
def is_bool(f):
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in ('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)'): return True
    return False
cand=set()
for a in CL: cand|=set(E.cone(a)[1])
cand-={18956,1530,1603}
cand=sorted(cand)
print(f"candidates {len(cand)} (non-boolean {sum(1 for f in cand if not is_bool(f))})",flush=True)
cols={}; knobs=[]
t0=time.time()
for f in cand:
    o=v0[f]
    r={}
    ok=True
    b1,_=fast.resid_delta(v0,bad0,{f:o+1})
    b2,_=fast.resid_delta(v0,bad0,{f:o+2})
    b7,_=fast.resid_delta(v0,bad0,{f:o+7})
    for a in CL:
        d1=b1.get(a,0)-bad0.get(a,0); d2=b2.get(a,0)-bad0.get(a,0); d7=b7.get(a,0)-bad0.get(a,0)
        if d2!=2*d1 or d7!=7*d1: ok=False; break
        if d1: r[a]=d1
    if ok and r:
        knobs.append(f); cols[f]=(r, set(b1)-set(bad0))
print(f"cluster-affine knobs: {len(knobs)} ({time.time()-t0:.0f}s)",flush=True)
for f in knobs:
    print("   x_%d %s: %s  (also disturbs %s)"%(f,'BOOL' if is_bool(f) else 'int',
          {a:str(c)[:20] for a,c in cols[f][0].items()}, sorted(cols[f][1])[:6]),flush=True)
rows=[{f:cols[f][0][a] for f in knobs if a in cols[f][0]} for a in CL]
rhs=[-bad0.get(a,0) for a in CL]
sol,msg,_=sparse.solve_sparse(rows,rhs,names=CL,verbose=True,maxcore=400)
print("CLUSTER SOLVE ->",msg,flush=True)
if sol is not None:
    ns=dict(s)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"EXACT: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)}",flush=True)
    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('close4_%d.json'%(39033-len(ff)),'w'))
    json.dump({str(k):str(int(x)) for k,x in ns.items()}, open('close4_seed.json','w'))
