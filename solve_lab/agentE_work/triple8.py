import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(5000000)
import engine as E, fast, harness as H, sparse
s={int(k):int(v) for k,v in json.load(open('triple_state_seed.json')).items()}
v0=E.forward(s); bad0=E.badatoms(v0)
def cf(a): return set(E.cone(a)[1])
cand=set()
for a in bad0: cand|=cf(a)
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
cand=sorted(cand-{18956,1530,1603})
cols={}; affine=[]
for f in cand:
    o=v0[f]
    b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2}); b7,_=fast.resid_delta(v0,bad0,{f:o+7})
    col={a:b1.get(a,0)-bad0.get(a,0) for a in set(b1)|set(bad0) if b1.get(a,0)-bad0.get(a,0)}
    if all(b2.get(a,0)-bad0.get(a,0)==2*col.get(a,0) and b7.get(a,0)-bad0.get(a,0)==7*col.get(a,0)
           for a in set(b2)|set(b7)|set(col)|set(bad0)):
        affine.append(f); cols[f]=col
atoms=sorted(set(bad0)|set().union(*[set(cols[f]) for f in affine]))
use=[a for a in atoms if a not in (20215,28647)]
rowmap={a:{} for a in use}
for f in affine:
    for a,c in cols[f].items():
        if a in rowmap: rowmap[a][f]=c
sol,msg,_=sparse.solve_sparse([rowmap[a] for a in use],[-bad0.get(a,0) for a in use],names=use,verbose=True,maxcore=400)
print("44-row affine solve ->",msg,flush=True)
if sol is not None:
    ns=dict(s)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"EXACT: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)}",flush=True)
    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('triple8_%d.json'%(39033-len(ff)),'w'))
    json.dump({str(k):str(int(x)) for k,x in ns.items()}, open('triple8_seed.json','w'))
