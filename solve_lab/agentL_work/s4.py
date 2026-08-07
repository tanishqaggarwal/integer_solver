"""Batched |S|=4 sweep: generic merge-by-deepest-LCA over 4 leaves."""
import importlib.util, itertools, time, pickle, sys
spec=importlib.util.spec_from_file_location('ss','/home/user/integer_solver/solve_lab/agentL_work/subsearch.py')
ss=importlib.util.module_from_spec(spec); spec.loader.exec_module(ss)
p=ss.p; K=ss.K; TGT=ss.TGT; live=ss.live; ORIENT=ss.ORIENT
depth=ss.depth; parent=ss.parent; swapup=ss.swapup; sw2root=ss.sw2root
LEAF=ss.LEAF; cums=ss.cums; path=ss.path
def batch_inv(xs):
    n=len(xs); pre=[1]*(n+1)
    for i,x in enumerate(xs): pre[i+1]=pre[i]*x%p
    inv=pow(pre[n],p-2,p); out=[0]*n
    for i in range(n-1,-1,-1):
        out[i]=inv*pre[i]%p; inv=inv*xs[i]%p
    return out
LC={}
for i,A in enumerate(live):
    for B in live[i+1:]:
        m=ss.lca(A,B); LC[(A,B)]=m; LC[(B,A)]=m
SWN={}
def swnode(a,b):
    k=(a,b)
    if k in SWN: return SWN[k]
    s=False; x=a
    while x!=b: s^=swapup[x]; x=parent[x]
    SWN[k]=s; return s
def sw(v,s): return (v[1],v[0]) if s else v
VA={}
for L in live:
    for i,n in enumerate(path[L]): VA[(L,n)]=sw(LEAF[L],cums[L][i])
# node-vs-node LCA via ancestor sets
anc={}
def ancs(n):
    if n in anc: return anc[n]
    s=set(); x=n
    while x!=ss.ROOT: s.add(x); x=parent[x]
    s.add(ss.ROOT); anc[n]=s; return s
def nlca(a,b):
    sa=ancs(a); x=b
    while x not in sa: x=parent[x]
    return x
pairs=list(itertools.combinations(live,2))
den=[]; meta=[]
for A,B in pairs:
    m=LC[(A,B)]; o=ORIENT[m]; X=VA[(A,m)]; Y=VA[(B,m)]
    den.append((Y[o]-X[o])%p or 1); meta.append((m,o,X,Y))
inv=batch_inv(den); PV={}
for (A,B),iv,(m,o,X,Y) in zip(pairs,inv,meta):
    ax,ay,bx,by=X[o],X[1-o],Y[o],Y[1-o]
    if (bx-ax)%p==0: PV[(A,B)]=PV[(B,A)]=None; continue
    l=(by-ay)*iv%p; ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
    PV[(A,B)]=PV[(B,A)]=(m,(ox,oy) if o==0 else (oy,ox))
print('pair table ready',flush=True)
CH=400000
def flush(buf,out):
    den=[(b[3][b[2]]-b[2+0][0] if False else 0) for b in []]
    d=[]; 
    for T,m,o,Aa,Bb in buf: d.append((Bb[o]-Aa[o])%p or 1)
    iv=batch_inv(d); res=[]
    for (T,m,o,Aa,Bb),v in zip(buf,iv):
        ax,ay,bx,by=Aa[o],Aa[1-o],Bb[o],Bb[1-o]
        if (bx-ax)%p==0: res.append((T,None,m)); continue
        l=(by-ay)*v%p; ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
        res.append((T,(ox,oy) if o==0 else (oy,ox),m))
    out.extend(res)
t0=time.time(); n=0; hits=[]
stage1=[]; stage2=[]
def process(chunk):
    """chunk: list of 4-tuples"""
    global hits
    # step A: merge the deepest-LCA pair of each quadruple using PV
    b1=[]
    for T in chunk:
        best=None
        for i in range(4):
            for j in range(i+1,4):
                m=LC[(T[i],T[j])]
                if best is None or depth[m]>best[0]: best=(depth[m],i,j,m)
        _,i,j,m1=best
        r=PV[(T[i],T[j])]
        if r is None: continue
        rest=[T[k] for k in range(4) if k not in (i,j)]
        # merge with whichever remaining leaf has the deeper LCA with m1
        c=[(depth[nlca(m1,z)],z) for z in rest]
        c.sort(reverse=True)
        _,Z=c[0]; W=[z for z in rest if z!=Z][0]
        m2=nlca(m1,Z)
        b1.append((T,m2,ORIENT[m2],sw(r[1],swnode(m1,m2)),VA[(Z,m2)],W))
    d=[(x[4][x[2]]-x[3][x[2]])%p or 1 for x in b1]
    iv=batch_inv(d); b2=[]
    for (T,m2,o,Aa,Bb,W),v in zip(b1,iv):
        ax,ay,bx,by=Aa[o],Aa[1-o],Bb[o],Bb[1-o]
        if (bx-ax)%p==0: continue
        l=(by-ay)*v%p; ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
        v2=(ox,oy) if o==0 else (oy,ox)
        m3=nlca(m2,W)
        b2.append((T,m3,ORIENT[m3],sw(v2,swnode(m2,m3)),VA[(W,m3)]))
    d=[(x[4][x[2]]-x[3][x[2]])%p or 1 for x in b2]
    iv=batch_inv(d)
    for (T,m3,o,Aa,Bb),v in zip(b2,iv):
        ax,ay,bx,by=Aa[o],Aa[1-o],Bb[o],Bb[1-o]
        if (bx-ax)%p==0: continue
        l=(by-ay)*v%p; ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
        v3=(ox,oy) if o==0 else (oy,ox)
        if sw(v3,sw2root[m3])==TGT: hits.append(T); print('HIT4',T,flush=True)
buf=[]
for T in itertools.combinations(live,4):
    buf.append(T); n+=1
    if len(buf)>=CH:
        process(buf); buf=[]
        print(n,'%.0fs'%(time.time()-t0),flush=True)
if buf: process(buf)
print('|S|=4 done: %d sets, %d hits, %.0fs'%(n,len(hits),time.time()-t0))
pickle.dump(hits,open('hits4.pkl','wb'))
