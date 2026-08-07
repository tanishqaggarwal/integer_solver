"""Exact minimum-violation (coset decoding) over a chosen departure support.

Every equation restricted to the support is an exact polynomial in |S| unknowns.  At the
base point (the all-linear-forced solution x0) the affine ones vanish and 20 constants do
not.  We look for the assignment minimising the number of nonzero equations, by testing,
for increasing k, whether some k-subset can be violated with the rest simultaneously
satisfiable.
"""
import os, sys, pickle, itertools, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']; n=len(NB); pt=D['pt']
ix={u:i for i,u in enumerate(NB)}
SUP=[ix[int(u)] for u in sys.argv[1].split(',')]
MAXK=int(sys.argv[2]) if len(sys.argv)>2 else 7
k=len(SUP); sidx={c:j for j,c in enumerate(SUP)}
print('support (%d unknowns): %s'%(k,[NB[c] for c in SUP]),flush=True)
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
polys=[]
for i,f in lin: polys.append((i,sub(f)))
for i,f in non: polys.append((i,sub(f if not isinstance(f,int) else {():f})))
varying=[(i,g) for i,g in polys if any(any(e) for e in g)]
constbad=[i for i,g in polys if g and not any(any(e) for e in g)]
print('varying %d ; nonzero-constant (unfixable) %d'%(len(varying),len(constbad)),flush=True)
if constbad: print('   unfixable:',constbad[:20])
aff=[(i,g) for i,g in varying if max(sum(m) for m in g)<=1]
hi=[(i,g) for i,g in varying if max(sum(m) for m in g)>1]
print('affine %d ; higher-degree %d %s'%(len(aff),len(hi),[i for i,_ in hi][:10]),flush=True)
base=[x0[c] for c in SUP]
def ev(g,p):
    s=0
    for m,c in g.items():
        t=c
        for j,e in enumerate(m):
            if e: t=t*pow(p[j],e,P)%P
        s=(s+t)%P
    return s
print('failing at base point: %d'%len([1 for i,g in varying if ev(g,base)]))
delta=[pt[c] for c in SUP]
print('failing at the deliverable point: %d -> %s'%(len([1 for i,g in varying if ev(g,delta)]),
      [i for i,g in varying if ev(g,delta)]),flush=True)
# --- exact search: which k-subsets of the affine equations can be violated ---
def rref(M,nc):
    M=[r[:] for r in M]; piv=[]; r=0
    for c in range(nc):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]%P: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        iv=pow(M[r][c],-1,P); M[r]=[x*iv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==len(M): break
    return M,piv,r
rows={}
for i,g in aff:
    row=[0]*(k+1)
    for m,c in g.items():
        if sum(m)==0: row[k]=(-c)%P
        else: row[[j for j,e in enumerate(m) if e][0]]=c%P
    rows[i]=row
def consistent(drop):
    M=[rows[i] for i in rows if i not in drop]
    MM,piv,r=rref(M,k)
    return not any(all(x%P==0 for x in MM[t][:k]) and MM[t][k]%P for t in range(len(MM)))
idxs=list(rows)
print('\naffine equations to satisfy: %d'%len(idxs),flush=True)
t0=time.time()
best=None
for kk in range(0,MAXK+1):
    hit=False
    for T in itertools.combinations(idxs,kk):
        if consistent(set(T)):
            print('  *** %d affine equations violated suffices: %s  (%.0fs)'%(kk,T,time.time()-t0),flush=True)
            best=T; hit=True; break
    if hit: break
    print('  no %d-subset works (%.0fs)'%(kk,time.time()-t0),flush=True)
print('minimum violated affine equations over this support:',len(best) if best is not None else '>%d'%MAXK)
