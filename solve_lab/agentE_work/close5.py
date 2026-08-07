"""Exact solve of the closed 5-row subsystem {20215,28647,7389,10187,20212} over the six
   integer knobs whose disturbance sets stay inside it."""
import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, sparse, harness as H
KN=[6083,11436,14393,14853,22820,26489,31339,37012]
s={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}
v0=E.forward(s); bad0=E.badatoms(v0)
print("start bad",sorted(bad0),flush=True)
cols={}; touched=set(bad0)
for f in KN:
    o=v0[f]
    b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2}); b7,_=fast.resid_delta(v0,bad0,{f:o+7})
    col={}
    aff=True
    for a in set(b1)|set(b2)|set(b7)|set(bad0):
        d1=b1.get(a,0)-bad0.get(a,0); d2=b2.get(a,0)-bad0.get(a,0); d7=b7.get(a,0)-bad0.get(a,0)
        if d2!=2*d1 or d7!=7*d1: aff=False
        if d1: col[a]=d1
    cols[f]=col; touched|=set(col)
    print(f"  x_{f}: affine={aff} touches {sorted(col)}",flush=True)
rowsA=sorted(touched)
print("rows:",rowsA,flush=True)
rows=[{f:cols[f][a] for f in KN if a in cols[f]} for a in rowsA]
rhs=[-bad0.get(a,0) for a in rowsA]
for a,r,b in zip(rowsA,rows,rhs):
    print(f"  a{a}: " + " + ".join(f"({str(c)[:14]}..)*d{f}" for f,c in r.items()) + f" = {str(b)[:24]}..",flush=True)
sol,msg,_=sparse.solve_sparse(rows,rhs,names=rowsA,verbose=True,maxcore=100)
print("5-ROW SUBSYSTEM ->",msg,flush=True)
if sol is not None:
    ns=dict(s)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"EXACT: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)}",flush=True)
    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('close6_%d.json'%(39033-len(ff)),'w'))
    json.dump({str(k):str(int(x)) for k,x in ns.items()}, open('close6_seed.json','w'))
