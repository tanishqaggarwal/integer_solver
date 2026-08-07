"""Restrict the higher-degree equations to the cost-free (affine-kernel) departure
directions and ask how many can be zeroed simultaneously."""
import os, sys, pickle, itertools, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
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
aff=[(i,g) for i,g in varying if max(sum(m) for m in g)<=1]
hi=[(i,g) for i,g in varying if max(sum(m) for m in g)>1]
rowsA=[]
for i,g in aff:
    r=[0]*k
    for m,c in g.items():
        if sum(m)==1: r[[j for j,e in enumerate(m) if e][0]]=c%P
    rowsA.append((i,r))
# kernel of the affine rows
M=[r[:] for _,r in rowsA]; piv=[]; rr=0
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
free=[c for c in range(k) if c not in piv]
print('affine rank %d ; cost-free kernel dimension %d ; kernel coordinates %s'%(rr,len(free),[NB[SUP[c]] for c in free]))
K=[]
for j,c in enumerate(free):
    vec=[0]*k; vec[c]=1
    for t,cp in enumerate(piv): vec[cp]=(-M[t][c])%P
    K.append(vec)
base=[x0[c] for c in SUP]
m=len(K)
def subK(g):
    """substitute d = base + sum t_j K[j]  -> polynomial in m unknowns"""
    out={}
    # expand each monomial
    for mono,c in g.items():
        cur={(0,)*m:c%P}
        for j,e in enumerate(mono):
            for _ in range(e):
                lin_j={(0,)*m:base[j]}
                for t in range(m):
                    ee=[0]*m; ee[t]=1
                    if K[t][j]%P: lin_j[tuple(ee)]=K[t][j]%P
                nxt={}
                for m1,c1 in cur.items():
                    for m2,c2 in lin_j.items():
                        mm=tuple(a+b for a,b in zip(m1,m2))
                        nxt[mm]=(nxt.get(mm,0)+c1*c2)%P
                cur={a:b for a,b in nxt.items() if b}
        for a,b in cur.items(): out[a]=(out.get(a,0)+b)%P
    return {a:b for a,b in out.items() if b}
HK=[(i,subK(g)) for i,g in hi]
nz=[(i,g) for i,g in HK if g]
print('higher-degree equations restricted to the kernel: %d nonzero, %d identically zero'%(len(nz),len(HK)-len(nz)))
for i,g in nz[:25]:
    print('   eq%-6d deg %d  %d terms'%(i,max(sum(a) for a in g),len(g)))
# span dimension of these polynomials
mons=sorted({a for _,g in nz for a in g})
mi={a:t for t,a in enumerate(mons)}
rows=[[g.get(a,0) for a in mons] for _,g in nz]
def rank(R):
    R=[r[:] for r in R]; rr=0
    for c in range(len(mons)):
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
print('monomials %d ; span dimension of the restricted higher-degree equations: %d'%(len(mons),rank(rows)))
pickle.dump({'K':K,'base':base,'HK':HK,'SUP':[NB[c] for c in SUP],'aff':rowsA},
            open('/home/user/integer_solver/solve_lab/agentG_work/kernelsys.pkl','wb'))
