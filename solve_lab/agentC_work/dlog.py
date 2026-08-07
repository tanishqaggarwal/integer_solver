"""Due diligence on the discrete log: small k, low Hamming weight, structured k."""
import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ec import *
chain=[int(b) for b in json.load(open('/home/user/integer_solver/solve_lab/agentC_work/chain.json'))]
pts=leafpoints()
G=pts[chain[0]]
C=json.load(open('/home/user/integer_solver/solve_lab/agentC_work/curve.json'))
Q=(int(C['Q'][0]),int(C['Q'][1]))
# sanity: chain[i] should be 2^i G
assert mul(2,pts[chain[0]])==pts[chain[1]]
print('chain[i]=2^i G verified for i=1')
print('chain[255]==2^255 G ?', mul(pow(2,255),G)==pts[chain[255]])
# 1) is Q a leaf?
S={pts[b]:i for i,b in enumerate(chain)}
print('Q is leaf?', S.get(Q))
# 2) BSGS in [0,2^44) with batch inversion
def batch_inv(xs):
    n=len(xs); pre=[1]*(n+1)
    for i,x in enumerate(xs): pre[i+1]=pre[i]*x%P
    inv=pow(pre[n],P-2,P)
    out=[0]*n
    for i in range(n-1,-1,-1):
        out[i]=pre[i]*inv%P; inv=inv*xs[i]%P
    return out
def table(base,m,start=None):
    """returns dict x-> j for j*base, j=0..m-1 (affine, batched)"""
    d={}; cur=start
    B=2048
    R=None
    j=0
    pts_=[]
    # simple incremental with batching
    acc=None
    xs=[]; buf=[]
    while j<m:
        k=min(B,m-j)
        # compute acc + i*base for i in 0..k-1 sequentially but batched inversion isn't trivial
        for i in range(k):
            if acc is None: d[0]=j
            else: d[acc[0]]=j
            acc=add(acc,base); j+=1
        if j%(1<<20)==0: print('   baby',j,time.time()-T0,flush=True)
    return d
T0=time.time()
M=1<<21
print('building baby table of size',M)
tb=table(G,M)
print('baby done',time.time()-T0,flush=True)
mg=mul(M,G); mgn=neg(mg)
R=Q; found=None
for i in range(M):
    if R is not None and R[0] in tb:
        j=tb[R[0]]
        for cand in (i*M+j, i*M-j):
            if mul(cand,G)==Q: found=cand; break
        if found: break
    R=add(R,mgn)
    if i%(1<<19)==0: print('   giant',i,time.time()-T0,flush=True)
print('BSGS range 2^42 result:',found)
json.dump({'small_k':found},open('/home/user/integer_solver/solve_lab/agentC_work/dlog_small.json','w'))
