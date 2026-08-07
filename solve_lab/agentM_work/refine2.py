"""Refinement, done right: two leaves stay together iff their class labels agree in every
   configuration in which NEITHER of them is ON."""
import sys,json,collections,random,time,pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentM_work')
import mcore as M, xcompare as X
s0=M.load_seed()
BASE=dict(s0); BASE[1530]=0; BASE[1603]=0
LEAVES=M.bools(); T,NODES=M.tree()

def labels(seed):
    v0,bad0,sig=M.measure(seed,coordfull=True)
    ids={}; out={}
    for f,s in sig.items():
        if s=='ON': out[f]=-1; continue
        if s not in ids: ids[s]=len(ids)
        out[f]=ids[s]
    return out,len(bad0),{v:k for k,v in ids.items()}

cfgs=[('alloff',dict(BASE))]
for L in LEAVES:
    c=dict(BASE); c[L]=1; cfgs.append(('on%d'%L,c))
# a few 2-leaf configs to break remaining ties
random.seed(5)
for r in range(40):
    c=dict(BASE)
    for L in random.sample(LEAVES,2): c[L]=1
    cfgs.append(('p%d'%r,c))

MAT=[]; TAGS=[]; SIGS=[]
t0=time.time()
for i,(tag,c) in enumerate(cfgs):
    lab,nb,inv=labels(c); MAT.append(lab); TAGS.append(tag); SIGS.append(inv)
    if i%40==0: print('%d/%d %s  t=%.0fs'%(i,len(cfgs),tag,time.time()-t0),flush=True)
pickle.dump({'MAT':MAT,'TAGS':TAGS,'SIGS':SIGS,'LEAVES':LEAVES},open('refine2.pkl','wb'))

# union-find on "compatible in every config where neither is ON"
par=list(range(len(LEAVES))); idx={f:i for i,f in enumerate(LEAVES)}
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
def uni(a,b):
    a,b=find(a),find(b)
    if a!=b: par[a]=b
for i in range(len(LEAVES)):
    for j in range(i+1,len(LEAVES)):
        u,v=LEAVES[i],LEAVES[j]
        if find(i)==find(j): continue
        ok=True
        for lab in MAT:
            lu,lv=lab[u],lab[v]
            if lu==-1 or lv==-1: continue
            if lu!=lv: ok=False; break
        if ok: uni(i,j)
blocks=collections.defaultdict(list)
for i,f in enumerate(LEAVES): blocks[find(i)].append(f)
B=sorted((sorted(v) for v in blocks.values()),key=len,reverse=True)
print('BLOCKS:',len(B),[len(b) for b in B])
pickle.dump(B,open('blocks2.pkl','wb'))
# compare each block to tree nodes
for b in B:
    s=set(b); nd,left=X.decompose(s)
    exact=[k for k,g in X.SUB if g==s]
    print(' n=%-4d exact_node=%s  cover=%s loose=%d'%(len(b),exact or '-',
          ','.join('%s(%d)'%(a,c) for a,c in nd) or '-',len(left)))
