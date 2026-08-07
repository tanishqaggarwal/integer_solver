"""Maximal common refinement of channel partitions over many configurations."""
import sys,json,collections,random,time,pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentM_work')
import mcore as M, xcompare as X

s0=M.load_seed()
BASE=dict(s0); BASE[1530]=0; BASE[1603]=0
LEAVES=M.bools()
T,NODES=M.tree(); SUB=X.SUB

def sigmap(seed):
    """leaf -> class label (full 5-coord signature); ON leaves get label None."""
    v0,bad0,sig=M.measure(seed,coordfull=True)
    return {f:(None if s=='ON' else s) for f,s in sig.items()}, len(bad0)

def run(cfgs,outfile):
    key={f:() for f in LEAVES}
    log=[]
    t0=time.time()
    for i,(tag,c) in enumerate(cfgs):
        sm,nb=sigmap(c)
        ncl=len(set(v for v in sm.values() if v is not None))
        for f in LEAVES:
            key[f]=key[f]+(sm[f],)
        log.append((tag,nb,ncl))
        if i%25==0: print('  %d/%d %s cls=%d  t=%.0fs'%(i,len(cfgs),tag,ncl,time.time()-t0),flush=True)
    blocks=collections.defaultdict(list)
    for f in LEAVES: blocks[key[f]].append(f)
    B=sorted((sorted(v) for v in blocks.values()),key=len,reverse=True)
    pickle.dump({'blocks':B,'log':log},open(outfile,'wb'))
    return B,log

if __name__=='__main__':
    cfgs=[('alloff',dict(BASE))]
    for L in LEAVES:
        c=dict(BASE); c[L]=1; cfgs.append(('on%d'%L,c))
    B,log=run(cfgs,'refine256.pkl')
    print('BLOCKS:',len(B),[len(b) for b in B])
