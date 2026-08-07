"""Fast forward-evaluation harness: free inputs -> all vars -> violated atoms."""
import pickle, re, collections, json, math, sys
m=pickle.load(open('/home/user/integer_solver/solve_lab/agentM_work/model3.pkl','rb'))
atoms=m['atoms']; eqt=m['eq_terms']
d=pickle.load(open('/home/user/integer_solver/solve_lab/agentM_work/dag.pkl','rb'))
info=d['info']; FREE=sorted(d['free'])
VAR=re.compile(r'x_(\d+)'); NV=38748
avars=[sorted(vs) for _,vs in info]
acodes=[compile(VAR.sub(r'v[\1]',a),'<a>','eval') for a in atoms]
occ=collections.defaultdict(list)
for i,vs in enumerate(avars):
    for x in vs: occ[x].append(i)

def _bootstrap():
    """One slow propagation from free=0 to fix the orientation (definer map + topo order)."""
    v=[0]*NV; known=[False]*NV; definer=[None]*NV
    for x in FREE: known[x]=True
    ns={'v':v,'__builtins__':{}}
    q=collections.deque(range(len(atoms))); inq=[True]*len(atoms); seq=[]
    while q:
        i=q.popleft(); inq[i]=False
        unk=[x for x in avars[i] if not known[x]]
        if len(unk)!=1: continue
        u=unk[0]; old=v[u]
        v[u]=0; c0=eval(acodes[i],ns); v[u]=1; c1=eval(acodes[i],ns); v[u]=2; c2=eval(acodes[i],ns); v[u]=old
        A2=c2-2*c1+c0
        if A2==0:
            sl=c1-c0
            if sl==0 or c0%sl: continue
            val=-c0//sl; kind=('lin',sl)
        else:
            A=A2//2; B=c1-c0-A; C=c0; disc=B*B-4*A*C
            if disc<0: continue
            r=math.isqrt(disc)
            if r*r!=disc: continue
            roots={(-B+s)//(2*A) for s in (r,-r) if (-B+s)%(2*A)==0}
            if len(roots)!=1: continue
            val=roots.pop(); kind=('quad',A)
        v[u]=val; known[u]=True; definer[u]=(i,kind); seq.append(u)
        for j in occ[u]:
            if not inq[j]: inq[j]=True; q.append(j)
    return definer, seq

import os
if os.path.exists('orient.pkl'):
    O=pickle.load(open('orient.pkl','rb')); definer=O['definer']; SEQ=O['seq']
else:
    definer,SEQ=_bootstrap(); pickle.dump({'definer':definer,'seq':SEQ}, open('orient.pkl','wb'))

# build explicit solved expression per var in SEQ, as a compiled closure
# var u defined by atom i: atom is linear or quadratic in u -> solve numerically each time.
SOLVE=[]
for u in SEQ:
    i,kind=definer[u]
    SOLVE.append((u,i,kind[0]))
CHECKS=[i for i in range(len(atoms)) if definer[u] is None or True]
defatom=set(i for u in SEQ for i,_ in [definer[u]])
CHECK_ATOMS=[i for i in range(len(atoms)) if i not in defatom]

def forward(freevals):
    v=[0]*NV
    for k,val in freevals.items(): v[k]=val
    ns={'v':v,'__builtins__':{}}
    for u,i,kd in SOLVE:
        c=acodes[i]
        v[u]=0; c0=eval(c,ns)
        v[u]=1; c1=eval(c,ns)
        if kd=='lin':
            sl=c1-c0
            v[u]=-c0//sl if sl and c0%sl==0 else 0
        else:
            v[u]=2; c2=eval(c,ns)
            A2=c2-2*c1+c0; A=A2//2; B=c1-c0-A; C=c0
            disc=B*B-4*A*C
            if disc<0: v[u]=0; continue
            r=math.isqrt(disc)
            if r*r!=disc: v[u]=0; continue
            rts={(-B+s)//(2*A) for s in (r,-r) if (-B+s)%(2*A)==0}
            v[u]=rts.pop() if len(rts)==1 else 0
    return v

def badatoms(v):
    ns={'v':v,'__builtins__':{}}
    return [i for i in range(len(atoms)) if eval(acodes[i],ns)!=0]

def eqfails(v):
    ns={'v':v,'__builtins__':{}}
    av={}
    bad=[]
    for i in range(len(atoms)):
        r=eval(acodes[i],ns)
        if r: av[i]=r
    fails=[]
    for e,(issq,outer,terms) in enumerate(eqt):
        s=0
        for c,a in terms:
            if a<0: s+=c
            elif a in av: s+=c*av[a]
        if s: fails.append(e)
    return fails,av

def dump(v,path):
    json.dump({f"x_{i}":v[i] for i in range(NV) if v[i]!=0}, open(path,'w'))
