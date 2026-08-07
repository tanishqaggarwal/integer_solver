"""Batched-inversion fold for |S|=2 and |S|=3, to make |S|=4 plausible.
Montgomery batch inversion: N inversions for the cost of 1 + 3N mults."""
import importlib.util, itertools, time, pickle
spec=importlib.util.spec_from_file_location('ss','/home/user/integer_solver/solve_lab/agentT_work/mirror/L/subsearch.py')
ss=importlib.util.module_from_spec(spec); spec.loader.exec_module(ss)
p=ss.p; K=ss.K; live=ss.live
def batch_inv(xs):
    n=len(xs); pre=[1]*(n+1)
    for i,x in enumerate(xs): pre[i+1]=pre[i]*x%p
    inv=pow(pre[n],p-2,p); out=[0]*n
    for i in range(n-1,-1,-1):
        out[i]=inv*pre[i]%p; inv=inv*xs[i]%p
    return out
t0=time.time()
xs=[(i*7919+13)%p or 1 for i in range(200000)]
iv=batch_inv(xs)
print('batch inversion of 200k: %.1fs  (correct: %s)'%(time.time()-t0, all(a*b%p==1 for a,b in list(zip(xs,iv))[:50])))
t0=time.time()
for x in xs[:20000]: pow(x,p-2,p)
print('naive pow inversion, per-inversion: %.0f us'%((time.time()-t0)/20000*1e6))
