"""Restricted to a departure support: the higher-degree equations' span dimension and
their proportionality classes.  If the span is 2-dimensional (a pencil in A and B), then
for any point with (A,B) != 0 the number that vanish equals the size of one class."""
import os, sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
W='/home/user/integer_solver/solve_lab/agentG_work/'
D=pickle.load(open(W+'coset_model.pkl','rb')); Lin=pickle.load(open(W+'coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']; pt=D['pt']
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
hi=[(i,g) for i,g in varying if max(sum(m) for m in g)>1]
base=[x0[c] for c in SUP]
def ev(g,p):
    s=0
    for m,c in g.items():
        t=c
        for j,e in enumerate(m):
            if e: t=t*pow(p[j],e,P)%P
        s=(s+t)%P
    return s
fail0=[i for i,g in hi if ev(g,base)]
print('support %d unknowns ; higher-degree equations %d ; nonzero at x0: %d'%(k,len(hi),len(fail0)))
mons=sorted({m for _,g in hi for m in g}); nm=len(mons)
V=[[g.get(m,0) for m in mons] for _,g in hi]
def rref(R):
    R=[r[:] for r in R]; piv=[]; rr=0
    for c in range(nm):
        pr=None
        for t in range(rr,len(R)):
            if R[t][c]%P: pr=t;break
        if pr is None: continue
        R[rr],R[pr]=R[pr],R[rr]
        iv=pow(R[rr][c],-1,P); R[rr]=[x*iv%P for x in R[rr]]
        for t in range(len(R)):
            if t!=rr and R[t][c]%P:
                f=R[t][c]; R[t]=[(x-f*y)%P for x,y in zip(R[t],R[rr])]
        piv.append(c); rr+=1
    return rr,piv,R
rr,piv,_=rref(V)
print('monomials %d ; SPAN DIMENSION of the higher-degree equations on this support: %d'%(nm,rr))
# proportionality classes
grp=collections.defaultdict(list)
for t,(i,g) in enumerate(hi):
    v=V[t]
    j0=next((c for c in range(nm) if v[c]%P),None)
    if j0 is None: grp[('ZERO',)].append(i); continue
    iv=pow(v[j0],-1,P)
    grp[tuple(x*iv%P for x in v)].append(i)
sizes=sorted((len(v) for v in grp.values()),reverse=True)
print('proportionality classes: %d ; sizes %s'%(len(grp),sizes))
mx=max(sizes)
print('LARGEST class: %d'%mx)
print('=> if the span is a 2-dim pencil, any point with the pencil value nonzero leaves')
print('   at least %d - %d = %d higher-degree equations failing.'%(len(fail0),mx,len(fail0)-mx))
for kk,v in sorted(grp.items(),key=lambda kv:-len(kv[1]))[:6]:
    print('   class of size %d : %s'%(len(v),v))
