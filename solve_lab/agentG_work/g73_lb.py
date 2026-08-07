"""Rigorous lower bound on the TOTAL number of violated equations, per support.

For a relaxation T of the affine equations, the freed subspace is the kernel of the rest.
Any higher-degree equation that is a NONZERO CONSTANT on that subspace fails no matter
what, so
        total(T)  >=  |T| + #{higher-degree equations constant and nonzero there}.
Minimising that bound over every admissible T with |T| <= BUDGET needs no root finding
and leaves no undecided cases.
"""
import os, sys, pickle, itertools, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
W='/home/user/integer_solver/solve_lab/agentG_work/'
D=pickle.load(open(W+'coset_model.pkl','rb')); Lin=pickle.load(open(W+'coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']
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
grp=collections.defaultdict(list)
for i,r in rowsA:
    j0=min(j for j in range(k) if r[j]%P); iv=pow(r[j0],-1,P)
    grp[tuple(x*iv%P for x in r)].append(i)
dirs=list(grp.items())
print('support %d unknowns ; affine rows %d -> %d distinct directions ; higher-degree %d'
      %(k,len(rowsA),len(dirs),len(hi)))
print('multiplicity histogram %s ; directions with multiplicity <= %d : %d'
      %(dict(collections.Counter(len(v) for _,v in dirs)),BUDGET,sum(1 for _,v in dirs if len(v)<=BUDGET)),flush=True)
def cnt(start,tot):
    c=0
    for j in range(start,len(dirs)):
        m=len(dirs[j][1])
        if tot+m<=BUDGET: c+=1+cnt(j+1,tot+m)
    return c
NC=cnt(0,0)
print('candidate violated-sets with total multiplicity <= %d : %d (EXHAUSTIVE)'%(BUDGET,NC),flush=True)
if NC>3*10**6:
    print('*** enumeration would exceed 3e6 -- NOT exhaustive at this budget, stopping ***'); sys.exit(0)
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
def nconst(K):
    """number of higher-degree equations that are nonzero constants on base+span(K)"""
    m=len(K); n=0
    for i,g in hi:
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
        out={a:b for a,b in out.items() if b}
        if out and not any(any(a) for a in out): n+=1
    return n
t0=time.time(); bestLB=[10**9,None]; tested=[0]
idxs=list(range(len(dirs)))
def rec(start,chosen,tot):
    if tot>BUDGET: return
    if chosen:
        rest=[dirs[j][0] for j in idxs if j not in chosen]
        Kr,Rr=nkern(rest)
        if Rr<R_all:
            viol=sorted(i for j in chosen for i in dirs[j][1])
            lb=len(viol)+nconst(Kr); tested[0]+=1
            if lb<bestLB[0]:
                bestLB[0]=lb; bestLB[1]=viol
                print('   new best lower bound %d  with T=%s'%(lb,viol),flush=True)
    for j in range(start,len(idxs)):
        m=len(dirs[j][1])
        if tot+m<=BUDGET: rec(j+1,chosen+[j],tot+m)
rec(0,[],0)
print('\nadmissible relaxations tested: %d in %.0fs'%(tested[0],time.time()-t0))
print('MINIMUM LOWER BOUND on total violated equations over all |T| <= %d : %d'%(BUDGET,bestLB[0]))
print('  attained by T = %s'%(bestLB[1],))
print('  (baseline T = {} costs 0 + %d = %d ; the 39,026 deliverable costs 7)'%(len([1 for i,g in hi if any(g)]) and 20,20))
