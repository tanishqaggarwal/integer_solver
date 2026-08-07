"""Maximal-sound per-boolean scan: for each boolean free input, flip it, symbolize
ALL non-boolean free inputs, solve the linear part, substitute -> exact residual."""
import os, sys, time, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G, gsolve
from gsym2 import L, ad, P
import gGclose
SRC=os.environ.get('SRC','/home/user/integer_solver/solve_lab/s10/AG_39013.json')
LO=int(sys.argv[1]); HI=int(sys.argv[2])
OUT=sys.argv[3] if len(sys.argv)>3 else 'bigscan_%d_%d.pkl'%(LO,HI)
v0=L.load(SRC); ad.fwd(v0,rounds=6)
FREE=[u for u in range(L.NVARS) if u not in L.definer]
BOOL=[u for u in FREE if gGclose.isbool(u)]
NB=[u for u in FREE if not gGclose.isbool(u)]
n=len(NB)
def analyse(v):
    val,skipped=G.build(v,NB,cap=6)
    rows=[];nzc=[]
    for a in G.check_atoms():
        f=G.evalatom(a,val,6)
        if isinstance(f,int):
            if f%P: nzc.append((a,f%P))
        else: rows.append((a,f))
    lin=[(a,f) for a,f in rows if G.deg(f)==1]
    non=[(a,f) for a,f in rows if G.deg(f)>1]
    sp=[]
    for a,f in lin:
        r={}
        for m,c in f.items():
            if not m: r[n]=(-c)%P
            else: r[m[0][0]]=c%P
        sp.append(r)
    piv,R=gsolve.sparse_rref(sp,n)
    inc=[i for i,r in enumerate(sp) if False]
    ninc=0; incchecks=[]
    for i,r in enumerate(R):
        cols=[c for c in r if c!=n]
        if not cols and r.get(n,0)%P: ninc+=1; incchecks.append(lin[i][0])
    need=set()
    for a,f in non:
        for m in f:
            for k,e in m: need.add(k)
    expr={}
    for k in need:
        if k in piv:
            r=R[piv[k]]; e={}
            cst=r.get(n,0)%P
            if cst: e[()]=cst
            for c,vv in r.items():
                if c!=n and c!=k: e[((c,1),)]=(-vv)%P
            expr[k]=e if e else 0
        else: expr[k]={((k,1),):1}
    res=[]
    for a,f in non:
        out=0
        for m,c in f.items():
            t=c%P
            for k,e in m:
                for _ in range(e):
                    t=G.pmul(t,expr[k])
                    if t==0: break
                if t==0: break
            if t!=0: out=G.padd(out,t)
        if not (isinstance(out,int) and out%P==0): res.append((a,out if isinstance(out,int) else 'POLY'))
    return dict(nzc=[a for a,_ in nzc],ninc=ninc,incchecks=incchecks,res=res,rank=len(piv),skipped=skipped)
base=analyse(list(v0))
print('BASE ninc=%d nzc=%d res=%s'%(base['ninc'],len(base['nzc']),[(a,str(g)[:10]) for a,g in base['res']]),flush=True)
out={}
t0=time.time()
for i,b in enumerate(BOOL[LO:HI]):
    w=list(v0); w[b]=1-w[b]; ad.fwd(w,rounds=8)
    try: r=analyse(w)
    except Exception as e:
        out[b]=('ERR',str(e)[:40]); continue
    out[b]=r
    if r['ninc']==0 and not r['nzc'] and not r['res']:
        print('*** FULL MOD-P SOLUTION at bit x%d'%b,flush=True)
    if (i%25)==0:
        print('  %d/%d %.0fs  (x%d: ninc=%d nzc=%d res=%d)'%(LO+i,HI,time.time()-t0,b,r['ninc'],len(r['nzc']),len(r['res'])),flush=True)
    if (i%100)==0: pickle.dump({'base':base,'out':out},open(OUT,'wb'))
pickle.dump({'base':base,'out':out},open(OUT,'wb'))
cnt=collections.Counter()
for b,r in out.items():
    if isinstance(r,tuple): cnt['ERR']+=1; continue
    cnt[(r['ninc'],len(r['nzc']),len(r['res']))]+=1
print('summary:',cnt.most_common(20))
