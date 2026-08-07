"""Propagation recording the definer of each var; then cone analysis of violated atoms."""
import pickle, re, collections, sys, json, time, math
m=pickle.load(open('model3.pkl','rb')); atoms=m['atoms']
d=pickle.load(open('dag.pkl','rb')); info=d['info']; free=d['free']
VAR=re.compile(r'x_(\d+)'); NV=38748
codes=[compile(VAR.sub(r'v[\1]',a),'<a>','eval') for a in atoms]
avars=[sorted(vs) for _,vs in info]
occ=collections.defaultdict(list)
for i,vs in enumerate(avars):
    for x in vs: occ[x].append(i)

def run(seed=None, order=None):
    v=[0]*NV; known=[False]*NV; definer=[None]*NV
    for x in free: known[x]=True
    if seed:
        for k,val in seed.items(): v[int(k[2:]) if isinstance(k,str) and k.startswith('x_') else int(k)]=int(val)
    ns={'v':v,'__builtins__':{}}
    def ev(i): return eval(codes[i],ns)
    q=collections.deque(order if order else range(len(atoms))); inq=[True]*len(atoms)
    seq=[]
    while q:
        i=q.popleft(); inq[i]=False
        unk=[x for x in avars[i] if not known[x]]
        if len(unk)!=1: continue
        u=unk[0]; old=v[u]
        v[u]=0; c0=ev(i); v[u]=1; c1=ev(i); v[u]=2; c2=ev(i); v[u]=old
        A2=c2-2*c1+c0
        if A2==0:
            sl=c1-c0
            if sl==0 or c0%sl: continue
            v[u]=-c0//sl
        else:
            A=A2//2; B=c1-c0-A; C=c0; disc=B*B-4*A*C
            if disc<0: continue
            r=math.isqrt(disc)
            if r*r!=disc: continue
            roots={ (-B+s)//(2*A) for s in (r,-r) if (-B+s)%(2*A)==0 }
            if len(roots)!=1: continue
            v[u]=roots.pop()
        known[u]=True; definer[u]=i; seq.append(u)
        for j in occ[u]:
            if not inq[j]: inq[j]=True; q.append(j)
    bad=[i for i in range(len(atoms)) if all(known[x] for x in avars[i]) and ev(i)!=0]
    undet=[i for i in range(len(atoms)) if not all(known[x] for x in avars[i])]
    return v,known,definer,seq,bad,undet

v,known,definer,seq,bad,undet=run()
print("bad",len(bad),bad,"undet",len(undet),"known",sum(known))
pickle.dump({'definer':definer,'seq':seq,'free':free},open('definer.pkl','wb'))
# cone of each bad atom back to free inputs
def cone(atomid):
    seen=set(); stack=list(avars[atomid]); res=set()
    while stack:
        u=stack.pop()
        if u in seen: continue
        seen.add(u)
        di=definer[u]
        if di is None: res.add(u); continue
        for w in avars[di]:
            if w!=u: stack.append(w)
    return seen,res
for b in bad:
    s,f=cone(b)
    print(f"atom {b}: {atoms[b][:80]} | cone {len(s)} vars, free inputs {len(f)}")
    print("   free:",sorted(f)[:60])
