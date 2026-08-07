"""EXHAUSTIVE minimum TOTAL violation on a departure support.

cost(T) = |T|  +  min over the subspace freed by dropping T of the number of
higher-degree equations that are nonzero.  A failing higher-degree equation costs +1,
it is not fatal -- that is the correction to g64/g66.
Baseline T = {} costs 0 + 20 = 20.  The 39,026 deliverable costs 7 + 0 = 7.
"""
import os, sys, pickle, itertools, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import flint
W='/home/user/integer_solver/solve_lab/agentG_work/'
D=pickle.load(open(W+'coset_model.pkl','rb')); Lin=pickle.load(open(W+'coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']; pt=D['pt']
ix={u:i for i,u in enumerate(NB)}
SUP=[ix[int(u)] for u in sys.argv[1].split(',')]
BUDGET=int(sys.argv[2]) if len(sys.argv)>2 else 6
k=len(SUP); sidx={c:j for j,c in enumerate(SUP)}
def sub(f):
    out={}
    for m,c in f.items():
        t=c%P; e=[0]*k
        for col,ee in m:
            if col in sidx: e[sidx[col]]+=ee
            else: t=t*pow(x0[col],ee,P)%P
        if t:
            key=tuple(e); out[key]=(out.get(key,0)+t)%P
    return {m:c for m,c in out.items() if c}
polys=[(i,sub(f)) for i,f in lin]+[(i,sub(f if not isinstance(f,int) else {():f})) for i,f in non]
varying=[(i,g) for i,g in polys if any(any(e) for e in g)]
aff=[(i,g) for i,g in varying if max(sum(m) for m in g)<=1]
hi=[(i,g) for i,g in varying if max(sum(m) for m in g)>1]
rowsA=[]
for i,g in aff:
    r=[0]*k
    for m,c in g.items():
        if sum(m)==1: r[[j for j,e in enumerate(m) if e][0]]=c%P
    rowsA.append((i,tuple(r)))
print('support %d unknowns ; affine rows %d ; higher-degree %d'%(k,len(rowsA),len(hi)),flush=True)
grp=collections.defaultdict(list)
for i,r in rowsA:
    j0=min(j for j in range(k) if r[j]%P); iv=pow(r[j0],-1,P)
    grp[tuple(x*iv%P for x in r)].append(i)
dirs=list(grp.items())
print('distinct directions %d ; multiplicity histogram %s'%(len(dirs),dict(collections.Counter(len(v) for _,v in dirs))),flush=True)
print('directions with multiplicity <= %d : %d'%(BUDGET,sum(1 for _,v in dirs if len(v)<=BUDGET)),flush=True)
def nkern(rows):
    M=[list(r) for r in rows]; piv=[]; rr=0
    for c in range(k):
        pr=None
        for t in range(rr,len(M)):
            if M[t][c]%P: pr=t;break
        if pr is None: continue
        M[rr],M[pr]=M[pr],M[rr]
        iv=pow(M[rr][c],-1,P); M[rr]=[x*iv%P for x in M[rr]]
        for t in range(len(M)):
            if t!=rr and M[t][c]%P:
                f=M[t][c]; M[t]=[(x-f*y)%P for x,y in zip(M[t],M[rr])]
        piv.append(c); rr+=1
    fr=[c for c in range(k) if c not in piv]
    K=[]
    for c in fr:
        vec=[0]*k; vec[c]=1
        for t,cp in enumerate(piv): vec[cp]=(-M[t][c])%P
        K.append(vec)
    return K,rr
K_all,R_all=nkern([r for _,r in rowsA])
base=[x0[c] for c in SUP]
ctx=flint.fmpz_mod_poly_ctx(P)
def minfail(K):
    """min number of nonzero higher-degree equations over base + span(K)."""
    m=len(K)
    def expand(g):
        out=collections.defaultdict(int)
        for mono,c in g.items():
            cur={(0,)*m:c%P}
            for j,e in enumerate(mono):
                for _ in range(e):
                    linj={(0,)*m:base[j]}
                    for t in range(m):
                        ee=[0]*m; ee[t]=1
                        if K[t][j]%P: linj[tuple(ee)]=K[t][j]%P
                    nxt=collections.defaultdict(int)
                    for m1,c1 in cur.items():
                        for m2,c2 in linj.items():
                            kk=tuple(a+b for a,b in zip(m1,m2)); nxt[kk]=(nxt[kk]+c1*c2)%P
                    cur={a:b for a,b in nxt.items() if b}
            for a,b in cur.items(): out[a]=(out[a]+b)%P
        return {a:b for a,b in out.items() if b}
    E=[(i,expand(g)) for i,g in hi]
    nz=[(i,e) for i,e in E if e]
    if not nz: return 0,[]
    constbad=[i for i,e in nz if not any(any(a) for a in e)]
    varb=[(i,e) for i,e in nz if any(any(a) for a in e)]
    if not varb: return len(constbad),constbad
    eff=[t for t in range(m) if any(a[t] for _,e in varb for a in e)]
    if len(eff)==1:
        t=eff[0]
        cands={0}
        for i,e in varb:
            deg=max(a[t] for a in e)
            co=[sum(c for a,c in e.items() if a[t]==d)%P for d in range(deg+1)]
            while co and co[-1]==0: co.pop()
            if len(co)>1:
                for r,_ in (ctx(co).roots() or []): cands.add(int(r))
        bestn=len(varb)+len(constbad); bestset=None
        for s in cands:
            f=[i for i,e in varb if sum(c*pow(s,a[t],P) for a,c in e.items())%P]
            if len(f)+len(constbad)<bestn: bestn=len(f)+len(constbad); bestset=constbad+f
        return bestn,bestset
    return None,constbad
t0=time.time(); best=(20,'baseline T={}')
print('baseline (T = {}): total 0 + %s'%str(minfail(K_all)[0]),flush=True)
results=[]
def rec(start,chosen,tot):
    if tot>BUDGET: return
    if chosen:
        rest=[dirs[j][0] for j in idxs if j not in chosen]
        Kr,Rr=nkern(rest)
        if Rr<R_all:
            viol=sorted(i for j in chosen for i in dirs[j][1])
            mf,fs=minfail(Kr)
            if mf is None:
                results.append((None,viol))
            else:
                tot2=len(viol)+mf
                if tot2<=BUDGET+2:
                    print('   T=%s (|T|=%d) + %d failing higher-degree = TOTAL %d %s'
                          %(viol,len(viol),mf,tot2,fs if fs else ''),flush=True)
                results.append((tot2,viol))
    for j in range(start,len(idxs)):
        m=len(dirs[j][1])
        if tot+m<=BUDGET: rec(j+1,chosen+[j],tot+m)
idxs=list(range(len(dirs)))
rec(0,[],0)
ok=[r for r in results if r[0] is not None]
und=[r for r in results if r[0] is None]
print('\nadmissible relaxations tested: %d (decided %d, undecided %d) in %.0fs'%(len(results),len(ok),len(und),time.time()-t0))
if ok:
    b=min(ok); print('BEST TOTAL over this support: %d  with T = %s'%(b[0],b[1]))
if und: print('undecided (>=3 effective parameters): %d, e.g. %s'%(len(und),und[0][1]))
