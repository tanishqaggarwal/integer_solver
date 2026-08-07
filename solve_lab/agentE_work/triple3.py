"""Solve the obstruction triple using ONLY the non-boolean integer handles on (U,V)."""
import sys, json, math, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(2000000)
import engine as E, fast, sparse, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
U_IDX,V_IDX=29210,8736
NB=[3401,4012,25710,28954,30468,33169,33177,34801,37856]   # non-boolean handles
def model(s, knobs):
    v0=E.forward(s); U0,V0=v0[U_IDX],v0[V_IDX]
    cols={}; nl=[]
    for f in knobs:
        v1,_=fast.apply_delta(v0,{f:v0[f]+1})
        v2,_=fast.apply_delta(v0,{f:v0[f]+2})
        a1=v1[U_IDX]-U0; b1=v1[V_IDX]-V0
        a2=v2[U_IDX]-U0; b2=v2[V_IDX]-V0
        if a2!=2*a1 or b2!=2*b1: nl.append(f)
        cols[f]=(a1,b1)
    return v0,U0,V0,cols,nl
def solve_once(s, knobs, verbose=True):
    v0,U0,V0,cols,nl=model(s,knobs)
    if verbose:
        print(f"  U0%p={U0%P}\n  V0%p={V0%P}\n  nonlinear knobs: {nl}")
        for f in knobs: print(f"    x_{f}: dU={cols[f][0]} dV={cols[f][1]}")
    W=-1
    r1={f:5002401*cols[f][0]+15322661*cols[f][1] for f in knobs if cols[f]!=(0,0)}
    b1=-(5002401*U0+15322661*V0)
    r2={f:cols[f][0] for f in knobs if cols[f][0]}; r2[W]=-P
    b2=-U0
    sol,msg,_=sparse.solve_sparse([r1,r2],[b1,b2],names=['exact','modp'],verbose=verbose,maxcore=200)
    return v0,sol,msg
if __name__=='__main__':
    s={int(k):int(v) for k,v in json.load(open('triple_state_seed.json')).items()}
    for it in range(8):
        print(f"--- iteration {it} ---",flush=True)
        v0,sol,msg=solve_once(s,NB,verbose=(it==0))
        print("  solve ->",msg,flush=True)
        if sol is None: break
        ns=dict(s)
        for f,d in sol.items():
            if f==-1 or not d: continue
            ns[f]=v0[f]+d
        v=E.forward(ns)
        U,V=v[U_IDX],v[V_IDX]
        ok1=(U%P==0); ok2=(5002401*U+15322661*V==0)
        av=E.badatoms(v); ff=E.eqfails(av)
        print(f"  after: p|U={ok1}  exactrow=0 -> {ok2}   fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:14]}",flush=True)
        s=ns
        json.dump({str(k):str(int(x)) for k,x in s.items()}, open('triple3_seed.json','w'))
        if ok1 and ok2: 
            json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open(f'triple3_{39033-len(ff)}.json','w'))
            print("  TRIPLE CONDITIONS MET",flush=True)
            break
