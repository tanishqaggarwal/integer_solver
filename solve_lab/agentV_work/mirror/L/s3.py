"""|S|=3 sweep with batched inversion."""
import pickle, sys, time, itertools
import importlib.util
spec=importlib.util.spec_from_file_location('ss','/home/user/integer_solver/solve_lab/agentT_work/mirror/L/subsearch.py')
ss=importlib.util.module_from_spec(spec); spec.loader.exec_module(ss)
p=ss.p; K=ss.K; TGT=ss.TGT; live=ss.live; ORIENT=ss.ORIENT
lca=ss.lca; val_at=ss.val_at; sw=ss.sw; sw2root=ss.sw2root; chord=ss.chord
path=ss.path; ancset=ss.ancset; posn=ss.posn; depth=ss.depth; LEAF=ss.LEAF; cums=ss.cums
parent=ss.parent
# precompute pairwise LCA for all pairs (256^2)
LC={}
for i,A in enumerate(live):
    for B in live[i+1:]:
        m=lca(A,B); LC[(A,B)]=m; LC[(B,A)]=m
def valnode_at(v,n_from,n_to):
    """value v given in frame n_from, express in frame n_to (n_to ancestor of n_from)"""
    s=False; x=n_from
    while x!=n_to: s^=ss.swapup[x]; x=parent[x]
    return sw(v,s)
def fold3(A,B,C):
    # find the pair with the deepest LCA
    ms=[(depth[LC[(A,B)]],A,B,C),(depth[LC[(A,C)]],A,C,B),(depth[LC[(B,C)]],B,C,A)]
    ms.sort(reverse=True)
    _,X,Y,Z=ms[0]
    m1=LC[(X,Y)]
    v1=chord(val_at(X,m1),val_at(Y,m1),ORIENT[m1])
    if v1 is None: return None
    m2=LC[(X,Z)]
    if depth[m2]>depth[LC[(Y,Z)]]: m2=LC[(Y,Z)]
    a=valnode_at(v1,m1,m2); b=val_at(Z,m2)
    v2=chord(a,b,ORIENT[m2])
    return None if v2 is None else sw(v2,sw2root[m2])
if __name__=='__main__':
    # validate against full fold
    import random
    spec2=importlib.util.spec_from_file_location('ff','/home/user/integer_solver/solve_lab/agentT_work/mirror/L/fastfold.py')
    ff=importlib.util.module_from_spec(spec2); spec2.loader.exec_module(ff)
    rnd=random.Random(4); ok=0;bad=0
    for _ in range(200):
        S=rnd.sample(live,3)
        if ff.fold(S)==fold3(*S): ok+=1
        else: bad+=1; print('MISMATCH',S)
    print('fold3 validation: match %d mismatch %d'%(ok,bad),flush=True)
    if bad: sys.exit(1)
    t0=time.time(); n=0; hits=[]
    for T in itertools.combinations(live,3):
        n+=1
        if fold3(*T)==TGT: hits.append(T); print('HIT3',T,flush=True)
        if n%200000==0: print(n,'%.0fs'%(time.time()-t0),flush=True)
    print('|S|=3 done: %d triples, %d hits, %.0fs'%(n,len(hits),time.time()-t0))
    pickle.dump(hits,open('hits3.pkl','wb'))
