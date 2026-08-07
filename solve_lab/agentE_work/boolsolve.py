"""For each boolean selector flip, re-propagate exactly, then run the 8-knob affine cluster
   solve on top.  Mixed 0/1 + integer: the flip is a real re-propagation, not a derivative."""
import sys, json, math, re, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, sparse, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
KN=[6083,11436,14393,14853,22820,26489,31339,37012]
res=pickle.load(open('boolknob.pkl','rb'))
base={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}
bits=[f for f,(d,c,cur) in res.items() if d[10187]%P or d[20212]%P]
print(f"testing {len(bits)} boolean flips",flush=True)
best=None
for n,f in enumerate(sorted(bits)):
    s=dict(base); s[f]=1
    v0=E.forward(s); bad0=E.badatoms(v0)
    cols={}; touched=set(bad0)
    for k in KN:
        o=v0[k]
        b1,_=fast.resid_delta(v0,bad0,{k:o+1}); b2,_=fast.resid_delta(v0,bad0,{k:o+2})
        col={}
        for a in set(b1)|set(bad0):
            d=b1.get(a,0)-bad0.get(a,0)
            if d: col[a]=d
        cols[k]=col; touched|=set(col)
    rowsA=sorted(touched)
    rows=[{k:cols[k][a] for k in KN if a in cols[k]} for a in rowsA]
    rhs=[-bad0.get(a,0) for a in rowsA]
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=rowsA,verbose=False,maxcore=120,maxcorebits=5_000_000)
    tag=''
    if sol is not None:
        ns=dict(s)
        for k,d in sol.items():
            if d: ns[k]=v0[k]+d
        v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
        tag=f" -> EXACT fails={len(ff)} score={39033-len(ff)} bad={sorted(av)}"
        if best is None or len(ff)<best[0]:
            best=(len(ff),dict(ns),sorted(av),f)
            json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('boolsolve_%d.json'%(39033-len(ff)),'w'))
            json.dump({str(a):str(int(b)) for a,b in ns.items()}, open('boolsolve_seed.json','w'))
    print(f"x_{f}: start bad={sorted(bad0)} rows={len(rowsA)} -> {msg[:60]}{tag}",flush=True)
print("BEST:",best[0] if best else None, best[3] if best else None, best[2] if best else None,flush=True)
