import os, sys, json, pickle, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G, gsolve
from gsym2 import L, ad, P
TAG=sys.argv[1] if len(sys.argv)>1 else 'base'
d=pickle.load(open('allsym_%s.pkl'%TAG,'rb'))
SYMS=d['syms']; rows=d['rows']; nzc=d['nzc']
ix={u:i for i,u in enumerate(SYMS)}; n=len(SYMS)
lin=[(a,f) for a,f in rows if G.deg(f)==1]
non=[(a,f) for a,f in rows if G.deg(f)>1]
print('linear %d nonlinear %d nzc %d'%(len(lin),len(non),len(nzc)),flush=True)
sp=[]
for a,f in lin:
    r={}
    for m,c in f.items():
        if not m: r[n]=(-c)%P
        else: r[m[0][0]]=c%P
    sp.append(r)
t0=time.time()
piv,R=gsolve.sparse_rref(sp,n,verbose=True)
inc=gsolve.inconsistent(sp,n,piv,R)
print('rank %d free %d inconsistent %d  (%.0fs)'%(len(piv),n-len(piv),len(inc),time.time()-t0),flush=True)
# variables appearing in nonlinear polys
need=set()
for a,f in non:
    for m in f:
        for k,e in m: need.add(k)
print('vars in nonlinear part: %d'%len(need))
# expression of each needed var: if pivot -> const - sum free*coef ; else itself
expr={}
for k in need:
    if k in piv:
        r=R[piv[k]]
        e={}
        cst=r.get(n,0)%P
        if cst: e[()]=cst
        for c,v in r.items():
            if c!=n and c!=k: e[((c,1),)]=(-v)%P
        expr[k]=e if e else 0
    else:
        expr[k]={((k,1),):1}
res=[]
for a,f in non:
    out=0
    for m,c in f.items():
        t=c%P
        for k,e in m:
            for _ in range(e):
                t=G.pmul(t,expr[k]); 
                if t==0: break
            if t==0: break
        if t!=0: out=G.padd(out,t)
    res.append((a,out))
nzres=[(a,g) for a,g in res if not (isinstance(g,int) and g%P==0)]
print('nonlinear residual after substitution: %d nonzero of %d'%(len(nzres),len(res)))
for a,g in nzres[:20]:
    print('   a%-6d neq=%-3d %s'%(a,len(L.atom2eq.get(a,{})), ('CONST %d'%(g%P)) if isinstance(g,int) else 'deg %d terms %d'%(G.deg(g),len(g))))
pickle.dump({'piv':piv,'R':R,'res':res,'SYMS':SYMS,'n':n,'inc':inc},open('bigsolve_%s.pkl'%TAG,'wb'))
