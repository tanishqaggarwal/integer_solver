"""EQUATION-level exact solve: symbolic value of every equation's inner sum, then
linear solve + substitution.  This is strictly weaker (more permissive) than requiring
every check atom to vanish."""
import os, sys, pickle, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F, gsolve
import gsym2 as G
from gsym2 import L, ad, P
arg=sys.argv[1] if len(sys.argv)>1 else '-'
FL=[int(x) for x in arg.split(',') if x] if arg!='-' else []
v=list(F.v0)
for b in FL: v[b]=1-v[b]
ad.fwd(v,rounds=8)
NB=F.NB; n=len(NB)
t0=time.time()
val,sk=G.build(v,NB,cap=6)
allsym=[G.evalatom(a,val,6) for a in range(L.NA)]
bad=[]
for i,(m,sq,co) in enumerate(L.eq_atoms):
    s=0
    for a,c in co.items():
        t=allsym[a]
        if isinstance(t,int):
            if t%P: s=G.padd(s,(c*t)%P)
        else: s=G.padd(s,{mm:(c*cc)%P for mm,cc in t.items()})
    if not (isinstance(s,int) and s%P==0): bad.append((i,s))
print('non-trivial equations: %d  (%.0fs)'%(len(bad),time.time()-t0),flush=True)
lin=[(i,s) for i,s in bad if not isinstance(s,int) and G.deg(s)==1]
non=[(i,s) for i,s in bad if isinstance(s,int) or G.deg(s)>1]
print('linear %d nonlinear/const %d'%(len(lin),len(non)),flush=True)
sp=[]
for i,s in lin:
    r={}
    for mm,c in s.items():
        if not mm: r[n]=(-c)%P
        else: r[mm[0][0]]=c%P
    sp.append(r)
piv,R=gsolve.sparse_rref(sp,n)
inc=[lin[i][0] for i,r in enumerate(R) if not [c for c in r if c!=n] and r.get(n,0)%P]
print('linear rank %d free %d inconsistent %d (%.0fs)'%(len(piv),n-len(piv),len(inc),time.time()-t0),flush=True)
need=set()
for i,s in non:
    if isinstance(s,int): continue
    for mm in s:
        for k,e in mm: need.add(k)
expr={}
for k in need:
    if k in piv:
        r=R[piv[k]]; e={}
        cst=r.get(n,0)%P
        if cst: e[()]=cst
        for c,vv in r.items():
            if c!=n and c!=k: e[((c,1),)]=(-vv)%P
        expr[k]=e if e else 0
    else: expr[k]={((k,1),):1}
res=[]
for i,s in non:
    if isinstance(s,int):
        if s%P: res.append((i,s%P))
        continue
    out=0
    for mm,c in s.items():
        t=c%P
        for k,e in mm:
            for _ in range(e):
                t=G.pmul(t,expr[k])
                if t==0: break
            if t==0: break
        if t!=0: out=G.padd(out,t)
    if not (isinstance(out,int) and out%P==0): res.append((i,out))
print('EQUATION-LEVEL residual after linear solve: %d equations (%.0fs)'%(len(res),time.time()-t0))
cc=[(i,g%P) for i,g in res if isinstance(g,int)]
pp=[(i,g) for i,g in res if not isinstance(g,int)]
print('  nonzero constants: %d ; symbolic: %d'%(len(cc),len(pp)))
vars_=sorted({NB[k] for i,g in pp for mm in g for k,_ in mm})
print('  vars in residual:',vars_)
print('  sample:',[(i,('C' if isinstance(g,int) else 'deg%d/%dt'%(G.deg(g),len(g)))) for i,g in res[:20]])
pickle.dump({'res':res,'NB':NB,'flip':FL,'inc':inc},open('eqsolve_%s.pkl'%(arg.replace(',','_')),'wb'))
