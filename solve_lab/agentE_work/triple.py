"""Targeted solve of the size>=2 obstruction triple:
       5002401*U + 15322661*V = 0   over Z,     p | U
   using only the free variables in the cones of x_25848, x_17317, x_18682, x_28841."""
import sys, json, pickle, math, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, fast, sparse, bitfeas2 as B
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
U_IDX, V_IDX = 29210, 8736
def load_state():
    s={int(k):int(v) for k,v in json.load(open('triple_state_seed.json')).items()}
    return s
def cone_free(u):
    c=set(); st=[u]; seen=set()
    while st:
        w=st.pop()
        if w in seen: continue
        seen.add(w); dv=E.definer[w]
        if dv is None: c.add(w); continue
        for z in E.avars[dv[0]]:
            if z!=w: st.append(z)
    return c
def model(s, freevars, verbose=True):
    v0=E.forward(s)
    U0,V0=v0[U_IDX],v0[V_IDX]
    cols={}; nonlin=[]
    for f in freevars:
        v1,_=fast.apply_delta(v0,{f:v0[f]+1})
        v2,_=fast.apply_delta(v0,{f:v0[f]+2})
        a1=v1[U_IDX]-U0; b1=v1[V_IDX]-V0
        a2=v2[U_IDX]-U0; b2=v2[V_IDX]-V0
        if a2!=2*a1 or b2!=2*b1: nonlin.append(f)
        if a1 or b1: cols[f]=(a1,b1)
    if verbose:
        print(f"  U0 mod p = {U0%P}\n  V0 mod p = {V0%P}")
        print(f"  movers: {len(cols)} of {len(freevars)}   nonlinear: {len(nonlin)}")
    return v0,U0,V0,cols,nonlin
if __name__=='__main__':
    s=load_state()
    tot=set()
    for u in (25848,17317,18682,28841): tot|=cone_free(u)
    tot=sorted(tot)
    print(f"free vars in the four cones: {len(tot)}")
    v0,U0,V0,cols,nonlin=model(s,tot)
    movers=sorted(cols)
    for f in movers:
        a,b=cols[f]
        print(f"   x_{f}: dU={str(a)[:34]} dV={str(b)[:34]}")
    # system:  row1 = exact ;  row2 = U - p*w = 0
    W=-1
    r1={f:5002401*cols[f][0]+15322661*cols[f][1] for f in movers}
    b1=-(5002401*U0+15322661*V0)
    r2={f:cols[f][0] for f in movers}; r2[W]=-P
    b2=-U0
    rows=[r1,r2]; rhs=[b1,b2]
    sol,msg,_=sparse.solve_sparse(rows,rhs,names=['exact','modp'],verbose=True,maxcore=200)
    print("targeted 2-row solve (all movers) ->",msg)
    if sol is None:
        lin=[f for f in movers if f not in set(nonlin)]
        print("retrying with only the LINEAR movers:",lin)
        r1b={f:5002401*cols[f][0]+15322661*cols[f][1] for f in lin}
        r2b={f:cols[f][0] for f in lin}; r2b[W]=-P
        sol,msg,_=sparse.solve_sparse([r1b,r2b],[b1,b2],names=['exact','modp'],verbose=True,maxcore=200)
        print("linear-only solve ->",msg)
    if sol is not None:
        print("  moves:",{k:str(x)[:26] for k,x in sol.items() if x and k!=W})
        ns=dict(s)
        for f,d in sol.items():
            if f==W or not d: continue
            ns[f]=v0[f]+d
        v=E.forward(ns)
        U,V=v[U_IDX],v[V_IDX]
        print("  after move: U mod p =",U%P," V mod p =",V%P," exactrow =",5002401*U+15322661*V)
        av=E.badatoms(v); ff=E.eqfails(av)
        print(f"  EXACT: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:15]}")
        json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('triple_move.json','w'))
        json.dump({str(k):str(int(x)) for k,x in ns.items()}, open('triple_move_seed.json','w'))
