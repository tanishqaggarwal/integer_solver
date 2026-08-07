"""Single-move scan over ALL free inputs of frame B, starting from the 39,026 witness."""
import frameB as FB, json, time, math, sys, pickle
from frameB import Frame, State
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
B=len(base.fails)
print('base failing',B,'score',base.score(),flush=True)
best=(B,None)
t0=time.time(); tested=0
DELTAS=[1,-1,2,-2,p,-p,2*p,-2*p]
hits=[]
for i,u in enumerate(fr.free):
    cur=base.fv.get(u,0)
    for d in DELTAS:
        g=base.clone().set_free({u:cur+d}); tested+=1
        n=len(g.fails)
        if n<B:
            hits.append((n,u,d)); print('  IMPROVE failing=%d  x_%d += %s'%(n,u,'p' if d==p else d),flush=True)
            if n<best[0]: best=(n,g.clone())
    if i%1500==0: print('  %d/%d  %.0fs  tested=%d'%(i,len(fr.free),time.time()-t0,tested),flush=True)
print('single-move scan: %d moves tested, %.0fs, improvements=%s'%(tested,time.time()-t0,hits[:10]))
# exact atom-solving moves over the supports of the 7 nonzero atoms
def roots(st,a,X,deg=3):
    b0=st.fv.get(X,0); ys=[]; t=st.clone()
    for k in range(deg+1):
        t.set_free({X:b0+k}); ys.append(t.av[a])
    dd=[ys[:]]
    for k in range(deg):
        q=dd[-1]; dd.append([q[i+1]-q[i] for i in range(len(q)-1)])
    c=[dd[k][0] for k in range(deg+1)]
    if all(x==0 for x in c): return []
    o=max(k for k in range(deg+1) if c[k]!=0)
    if o==0: return []
    if o==1: return [] if c[0]%c[1] else [b0-c[0]//c[1]]
    if o==2:
        aa=c[2]; bb=2*c[1]-c[2]; cc=2*c[0]; disc=bb*bb-4*aa*cc
        if disc<0: return []
        r=math.isqrt(disc)
        if r*r!=disc: return []
        return [b0+s//(2*aa) for s in (-bb+r,-bb-r) if s%(2*aa)==0]
    return []
t0=time.time(); cnt=0
for a in base.nz():
    for X in fr.SUPV[a]:
        for rt in roots(base,a,X):
            if base.fv.get(X,0)==rt: continue
            g=base.clone().set_free({X:rt}); cnt+=1
            if len(g.fails)<B:
                hits.append((len(g.fails),X,'root')); print('  IMPROVE(root) failing=%d x_%d'%(len(g.fails),X),flush=True)
                if len(g.fails)<best[0]: best=(len(g.fails),g.clone())
print('exact-solve moves: %d tested %.0fs'%(cnt,time.time()-t0))
print('BEST failing',best[0])
if best[1] is not None and best[0]<B:
    st=best[1]
    out='H_%d_climb.json'%(39033-best[0])
    json.dump({('x_%d'%i):st.v[i] for i in range(38748) if st.v[i]!=0},open(out,'w'))
    print('WROTE',out)
