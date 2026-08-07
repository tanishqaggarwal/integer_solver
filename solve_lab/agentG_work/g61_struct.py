"""Structure of the departure problem on a support: affine row rank, and the
higher-degree equations that carry the 20 base failures."""
import os, sys, pickle, itertools, collections, random, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']; pt=D['pt']; n=len(NB)
ix={u:i for i,u in enumerate(NB)}
SUP=[ix[int(u)] for u in sys.argv[1].split(',')]
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
base=[x0[c] for c in SUP]; delta=[pt[c] for c in SUP]
def ev(g,p):
    s=0
    for m,c in g.items():
        t=c
        for j,e in enumerate(m):
            if e: t=t*pow(p[j],e,P)%P
        s=(s+t)%P
    return s
print('support %s'%[NB[c] for c in SUP])
print('affine %d, higher-degree %d'%(len(aff),len(hi)))
print('higher-degree equations (deg, #terms, value at base, value at deliverable):')
for i,g in hi:
    print('   eq%-6d deg %d  %2d terms  base=%s  deliv=%s'%(i,max(sum(m) for m in g),len(g),
        'NONZERO' if ev(g,base) else 'zero','NONZERO' if ev(g,delta) else 'zero'))
# affine rows: homogeneous part (their constant term is 0 at the base by construction)
rowsA=[]
for i,g in aff:
    r=[0]*k
    for m,c in g.items():
        if sum(m)==1: r[[j for j,e in enumerate(m) if e][0]]=c%P
    rowsA.append((i,r))
def rank(rows):
    M=[r[:] for r in rows]; rr=0
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
        rr+=1
    return rr
R=[r for _,r in rowsA]
print('\naffine homogeneous rank: %d of %d unknowns'%(rank(R),k))
# deliverable departure direction and how many affine rows it kills
dd=[(pt[c]-x0[c])%P for c in SUP]
kill=[i for (i,r) in rowsA if sum(a*b for a,b in zip(r,dd))%P==0]
print('deliverable direction kills %d of %d affine rows; violates %d'%(len(kill),len(rowsA),len(rowsA)-len(kill)))
print('violated:',[i for (i,r) in rowsA if sum(a*b for a,b in zip(r,dd))%P])
