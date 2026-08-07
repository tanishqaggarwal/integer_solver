import pickle, re, ast, collections, sys, json, time, math
m=pickle.load(open('model3.pkl','rb'))
atoms=m['atoms']
d=pickle.load(open('dag.pkl','rb')); info=d['info']
VAR=re.compile(r'x_(\d+)')
NV=38748
codes=[compile(VAR.sub(r'v[\1]',a),'<a>','eval') for a in atoms]
avars=[sorted(vs) for _,vs in info]
v=[0]*NV
known=[False]*NV
ns={'v':v,'__builtins__':{}}
def ev(i): return eval(codes[i],ns)
# occurrence index
occ=collections.defaultdict(list)
for i,vs in enumerate(avars):
    for x in vs: occ[x].append(i)
def unknown_count(i):
    return sum(1 for x in avars[i] if not known[x])
t0=time.time()
queue=collections.deque(range(len(atoms)))
inq=[True]*len(atoms)
violations=set()
solved=0
while queue:
    i=queue.popleft(); inq[i]=False
    unk=[x for x in avars[i] if not known[x]]
    if len(unk)==0:
        if ev(i)!=0: violations.add(i)
        else: violations.discard(i)
        continue
    if len(unk)!=1: continue
    u=unk[0]
    old=v[u]
    v[u]=0; c0=ev(i)
    v[u]=1; c1=ev(i)
    v[u]=2; c2=ev(i)
    v[u]=old
    a2=(c2-2*c1+c0)//2 if (c2-2*c1+c0)%2==0 else None
    if c2-2*c1+c0==0:
        # linear c0 + (c1-c0)*u
        sl=c1-c0
        if sl==0:
            if c0!=0: violations.add(i)
            continue
        if c0% sl!=0: continue   # non-integral -> skip (could still be resolved later)
        val=-c0//sl
        v[u]=val; known[u]=True; solved+=1
    else:
        # quadratic A u^2 + B u + C
        A=(c2-2*c1+c0)//2; B=c1-c0-A; C=c0
        disc=B*B-4*A*C
        if disc<0: violations.add(i); continue
        r=math.isqrt(disc)
        if r*r!=disc: continue
        roots=set()
        for s in (r,-r):
            if (-B+s)%(2*A)==0: roots.add((-B+s)//(2*A))
        if len(roots)!=1: continue
        v[u]=roots.pop(); known[u]=True; solved+=1
    for j in occ[u]:
        if not inq[j]: inq[j]=True; queue.append(j)
print(f"propagated: solved={solved} known={sum(known)} violations={len(violations)} t={time.time()-t0:.1f}s")
nzatoms=[i for i in range(len(atoms)) if unknown_count(i)==0 and ev(i)!=0]
print("determined-and-violated atoms:",len(nzatoms))
print("free/unknown vars:", NV-sum(known))
json.dump({f"x_{i}":v[i] for i in range(NV) if v[i]!=0}, open('prop_assign.json','w'))
pickle.dump({'v':v,'known':known}, open('prop.pkl','wb'))
