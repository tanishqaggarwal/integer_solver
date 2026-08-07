"""Sweep configurations; record channel partitions; test tree-compatibility; build common refinement."""
import sys,json,collections,random,time,pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentM_work')
import mcore as M, xcompare as X

s0=M.load_seed()
LEAVES=M.bools()
assert len(LEAVES)==256
T,NODES=M.tree()
SUB=X.SUB

def partition(seed):
    v0,bad0,sig=M.measure(seed,coordfull=True)
    cls=M.classes(sig)
    out=[]
    for k,v in cls.items(): out.append((k,frozenset(v)))
    return out,len(bad0)

def cross_count(part):
    n=0; det=[]
    for k,g in SUB:
        hit=[len(g&s) for _,s in part if g&s]
        if len(hit)>1: n+=1; det.append((k,len(g),sorted(hit,reverse=True)))
    return n,det

random.seed(11)
CFGS=[('cfg0',dict(s0))]
# off the two saturators
c=dict(s0); c[1530]=0; c[1603]=0; CFGS.append(('all_off',c))
# single leaf on (from each class + random)
for L in [16214,4701,490,6821,17760,24601,2081,47,91,438,1530,1603,35979,21074,12054,30448]:
    c=dict(s0); c[1530]=0; c[1603]=0; c[L]=1; CFGS.append(('one_%d'%L,c))
# random k-subsets on
for k in (2,3,5,8,20,60):
    for r in range(4):
        c=dict(s0); c[1530]=0; c[1603]=0
        for L in random.sample(LEAVES,k): c[L]=1
        CFGS.append(('rand%d_%d'%(k,r),c))

res=[]
t0=time.time()
for tag,c in CFGS:
    part,nb=partition(c)
    nc,det=cross_count(part)
    sizes=sorted((len(s) for _,s in part),reverse=True)
    lab=[('ON' if k=='ON' else 'IN' if k=='INERT' else str(len(s))) for k,s in sorted(part,key=lambda x:-len(x[1]))]
    res.append((tag,nb,sizes,lab,nc,det,[sorted(s) for _,s in part],[k if isinstance(k,str) else 'ch' for k,s in part]))
    print('%-12s bad0=%-3d parts=%s cross=%d %s'%(tag,nb,lab,nc,det[:4]))
print('elapsed',time.time()-t0)
pickle.dump(res,open('sweep1.pkl','wb'))
