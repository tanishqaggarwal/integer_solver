"""Sharp step: if every subset of the failing higher-degree equations of size
20 - j (j <= BUDGET-1) already spans the whole higher-degree space, then vanishing of
that many forces ALL of them to vanish.  Then #failing is 0 or >= BUDGET, and since
#failing = 0 needs |T| >= 7 (exhaustively established), total >= 7 everywhere."""
import os, sys, pickle, itertools, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
W='/home/user/integer_solver/solve_lab/agentG_work/'
D=pickle.load(open(W+'coset_model.pkl','rb')); Lin=pickle.load(open(W+'coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']
ix={u:i for i,u in enumerate(NB)}
SUP=[ix[int(u)] for u in sys.argv[1].split(',')]
JMAX=int(sys.argv[2]) if len(sys.argv)>2 else 5
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
F=[(i,g) for i,g in hi if ev(g,base)]
print('support %d unknowns ; higher-degree %d ; failing at x0: %d'%(k,len(hi),len(F)))
mons=sorted({m for _,g in F for m in g}); nm=len(mons)
V=[[g.get(m,0) for m in mons] for _,g in F]
def rank(rows):
    R=[r[:] for r in rows]; rr=0
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
        rr+=1
    return rr
FULL=rank(V)
print('span dimension of the %d failing higher-degree equations: %d (over %d monomials)'%(len(F),FULL,nm))
t0=time.time(); worst=None
for j in range(1,JMAX+1):
    mn=FULL
    for drop in itertools.combinations(range(len(F)),j):
        r=rank([V[t] for t in range(len(F)) if t not in drop])
        if r<mn: mn=r; worst=(j,drop,r)
    print('   removing any %d of them: minimum remaining span dimension = %d %s'
          %(j,mn,'(FULL - vanishing of the rest forces all to vanish)' if mn==FULL else '<-- DROPS'),flush=True)
    if mn<FULL: break
print('%.0fs'%(time.time()-t0))
if worst is None:
    print('\nCONCLUSION: any %d of the %d vanishing forces ALL of them to vanish.'%(len(F)-JMAX,len(F)))
    print('So the number failing is 0 or >= %d.  #failing = 0 needs |T| >= 7 (exhaustive),'%(JMAX+1))
    print('and #failing >= %d gives total >= |T| + %d >= %d.  Hence total >= 7 on this support.'%(JMAX+1,JMAX+1,JMAX+1+1))
else:
    print('\nspan drops when removing %s -> subset %s can vanish without the rest'%(worst[0],worst[1]))
