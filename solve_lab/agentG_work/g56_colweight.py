"""Coset decoding, step 2: how many equations a departure along each unknown disturbs,
and what the deliverable's own departure looks like in unknown-space."""
import os, sys, pickle, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsolve
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; pt=D['pt']; n=len(NB)
rows=[]
for i,f in lin:
    r={}
    for m,c in f.items():
        if m: r[m[0][0]]=c%P
        else: r[n]=(-c)%P
    rows.append((i,r))
col=collections.defaultdict(list)
for i,r in rows:
    for c in r:
        if c!=n: col[c].append(i)
w=sorted(((len(v),c) for c,v in col.items()))
print('unknowns that occur in at least one linear equation: %d of %d'%(len(col),n))
print('smallest per-unknown equation counts (an upper bound on the cost of moving it alone):')
for k,c in w[:25]: print('   x%-6d occurs in %d linear equations %s'%(NB[c],k,col[c][:8]))
hist=collections.Counter(k for k,_ in w)
print('histogram of per-unknown equation counts (first 15):',dict(sorted(hist.items())[:15]))
# the all-linear-forced solution x0 and the deliverable's departure
sp=[dict(r) for _,r in rows]
piv,R=gsolve.sparse_rref(sp,n)
print('\nlinear system: %d equations, rank %d, kernel dimension %d'%(len(sp),len(piv),n-len(piv)))
t={c:0 for c in range(n) if c not in piv}
x0=[0]*n
for c in range(n):
    if c in piv:
        r=R[piv[c]]; val=r.get(n,0)%P
        for c2,vv in r.items():
            if c2!=n and c2!=c: val=(val-vv*t.get(c2,0))%P
        x0[c]=val
    else: x0[c]=0
d=[(pt[c]-x0[c])%P for c in range(n)]
supp=[c for c in range(n) if d[c]]
print('deliverable departure d* = deliverable - x0 : support %d unknowns'%len(supp))
print('   (these are only meaningful modulo the kernel; the disturbed equations are what count)')
dist=[i for i,r in rows if sum(r.get(c,0)*d[c] for c in r if c!=n)%P]
print('equations disturbed by d*: %d -> %s'%(len(dist),dist))
pickle.dump({'colw':{NB[c]:len(v) for c,v in col.items()},'col':col,'rows':rows,'piv':piv,'R':R,'x0':x0},
            open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','wb'))
