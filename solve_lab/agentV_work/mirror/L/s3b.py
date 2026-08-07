"""Batched |S|=3 sweep (Montgomery batch inversion)."""
import importlib.util, itertools, time, pickle, sys
spec=importlib.util.spec_from_file_location('ss','/home/user/integer_solver/solve_lab/agentT_work/mirror/L/subsearch.py')
ss=importlib.util.module_from_spec(spec); spec.loader.exec_module(ss)
p=ss.p; K=ss.K; TGT=ss.TGT; live=ss.live; ORIENT=ss.ORIENT
depth=ss.depth; parent=ss.parent; swapup=ss.swapup; sw2root=ss.sw2root
LEAF=ss.LEAF; cums=ss.cums; posn=ss.posn; path=ss.path
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
    """swap flag from node a up to ancestor b"""
    k=(a,b)
    if k in SWN: return SWN[k]
    s=False; x=a
    while x!=b: s^=swapup[x]; x=parent[x]
    SWN[k]=s; return s
def sw(v,s): return (v[1],v[0]) if s else v
VA={}
for L in live:
    for i,n in enumerate(path[L]): VA[(L,n)]=sw(LEAF[L],cums[L][i])
# pass 1: all pair folds at their LCA
pairs=list(itertools.combinations(live,2))
t0=time.time()
den=[]; meta=[]
for A,B in pairs:
    m=LC[(A,B)]; o=ORIENT[m]; X=VA[(A,m)]; Y=VA[(B,m)]
    den.append((Y[o]-X[o])%p or 1); meta.append((m,o,X,Y))
inv=batch_inv(den)
PV={}
for (A,B),iv,(m,o,X,Y) in zip(pairs,inv,meta):
    ax,ay,bx,by=X[o],X[1-o],Y[o],Y[1-o]
    if (bx-ax)%p==0: PV[(A,B)]=None; continue
    l=(by-ay)*iv%p; ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
    v=(ox,oy) if o==0 else (oy,ox)
    PV[(A,B)]=(m,v); PV[(B,A)]=(m,v)
print('pair folds: %d in %.1fs'%(len(pairs),time.time()-t0),flush=True)
hits=[h for (A,B),r in PV.items() if r and sw(r[1],sw2root[r[0]])==TGT]
print('|S|=2 hits:',len(hits))
# pass 2: triples, chunked
t0=time.time(); CH=400000; buf=[]; n=0; hits3=[]
def flush(buf):
    den=[]; ok=[]
    for T,m2,o,Aa,Bb in buf:
        d=(Bb[o]-Aa[o])%p
        den.append(d or 1); ok.append(d!=0)
    inv=batch_inv(den)
    for (T,m2,o,Aa,Bb),iv,g in zip(buf,inv,ok):
        if not g: continue
        ax,ay,bx,by=Aa[o],Aa[1-o],Bb[o],Bb[1-o]
        l=(by-ay)*iv%p; ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
        v=(ox,oy) if o==0 else (oy,ox)
        if sw(v,sw2root[m2])==TGT: hits3.append(T); print('HIT3',T,flush=True)
for T in itertools.combinations(live,3):
    A,B,C=T
    dab,dac,dbc=depth[LC[(A,B)]],depth[LC[(A,C)]],depth[LC[(B,C)]]
    if dab>=dac and dab>=dbc: X,Y,Z=A,B,C
    elif dac>=dbc: X,Y,Z=A,C,B
    else: X,Y,Z=B,C,A
    r=PV[(X,Y)]
    if r is None: continue
    m1,v1=r
    m2=LC[(X,Z)]
    if depth[LC[(Y,Z)]]<depth[m2]: m2=LC[(Y,Z)]
    buf.append((T,m2,ORIENT[m2],sw(v1,swnode(m1,m2)),VA[(Z,m2)]))
    n+=1
    if len(buf)>=CH:
        flush(buf); buf=[]
        print(n,'%.0fs'%(time.time()-t0),flush=True)
if buf: flush(buf)
print('|S|=3 done: %d triples, %d hits, %.0fs'%(n,len(hits3),time.time()-t0))
pickle.dump(hits3,open('hits3.pkl','wb'))
