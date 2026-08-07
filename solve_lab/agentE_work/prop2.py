import pickle, re, collections, sys, json, time, math
m=pickle.load(open('model3.pkl','rb')); atoms=m['atoms']
d=pickle.load(open('dag.pkl','rb')); info=d['info']; free=d['free']
VAR=re.compile(r'x_(\d+)'); NV=38748
codes=[compile(VAR.sub(r'v[\1]',a),'<a>','eval') for a in atoms]
avars=[sorted(vs) for _,vs in info]
seed=json.load(open(sys.argv[1])) if len(sys.argv)>1 else {}
v=[0]*NV; known=[False]*NV
for x in free: known[x]=True
for k,val in seed.items():
    i=int(k[2:]) if k.startswith('x_') else int(k); v[i]=int(val)
ns={'v':v,'__builtins__':{}}
def ev(i): return eval(codes[i],ns)
occ=collections.defaultdict(list)
for i,vs in enumerate(avars):
    for x in vs: occ[x].append(i)
t0=time.time()
queue=collections.deque(range(len(atoms))); inq=[True]*len(atoms)
solved=0; noninteg=0
while queue:
    i=queue.popleft(); inq[i]=False
    unk=[x for x in avars[i] if not known[x]]
    if len(unk)!=1: continue
    u=unk[0]; old=v[u]
    v[u]=0; c0=ev(i); v[u]=1; c1=ev(i); v[u]=2; c2=ev(i); v[u]=old
    A2=c2-2*c1+c0
    if A2==0:
        sl=c1-c0
        if sl==0: continue
        if c0%sl!=0: noninteg+=1; continue
        v[u]=-c0//sl
    else:
        A=A2//2; B=c1-c0-A; C=c0
        disc=B*B-4*A*C
        if disc<0: continue
        r=math.isqrt(disc)
        if r*r!=disc: continue
        roots=set()
        for s in (r,-r):
            if (-B+s)%(2*A)==0: roots.add((-B+s)//(2*A))
        if len(roots)!=1: continue
        v[u]=roots.pop()
    known[u]=True; solved+=1
    for j in occ[u]:
        if not inq[j]: inq[j]=True; queue.append(j)
print(f"solved={solved} known={sum(known)}/{NV} noninteg={noninteg} t={time.time()-t0:.1f}s")
det=[i for i in range(len(atoms)) if all(known[x] for x in avars[i])]
bad=[i for i in det if ev(i)!=0]
print(f"fully-determined atoms {len(det)}, violated {len(bad)}")
json.dump({f"x_{i}":v[i] for i in range(NV) if v[i]!=0}, open(sys.argv[2] if len(sys.argv)>2 else 'prop2.json','w'))
pickle.dump({'v':v,'known':known,'bad':bad},open('prop2.pkl','wb'))
for i in bad[:30]: print("  BAD",i,atoms[i][:120],"=",str(ev(i))[:60])
