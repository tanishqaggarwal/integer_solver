"""Full exact reduce for arbitrary boolean flip sets (maximal sound symbol set)."""
import os, sys, time, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G, gsolve
from gsym2 import L, ad, P
import gGclose
SRC=os.environ.get('SRC','/home/user/integer_solver/solve_lab/s10/AG_39013.json')
v0=L.load(SRC); ad.fwd(v0,rounds=6)
FREE=[u for u in range(L.NVARS) if u not in L.definer]
NB=[u for u in FREE if not gGclose.isbool(u)]
n=len(NB)
def analyse(v, verbose=False):
    val,skipped=G.build(v,NB,cap=6)
    rows=[];nzc=[]
    for a in G.check_atoms():
        f=G.evalatom(a,val,6)
        if isinstance(f,int):
            if f%P: nzc.append((a,f%P))
        else: rows.append((a,f))
    lin=[(a,f) for a,f in rows if G.deg(f)==1]
    non=[(a,f) for a,f in rows if G.deg(f)>1]
    m=len(lin)
    sp=[]
    for a,f in lin:
        r={}
        for mm,c in f.items():
            if not mm: r[n]=(-c)%P
            else: r[mm[0][0]]=c%P
        sp.append(r)
    piv,R=gsolve.sparse_rref(sp,n)
    incchecks=[]
    for i,r in enumerate(R):
        cols=[c for c in r if c!=n]
        if not cols and r.get(n,0)%P: incchecks.append(lin[i][0])
    need=set()
    for a,f in non:
        for mm in f:
            for k,e in mm: need.add(k)
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
        for mm,c in f.items():
            t=c%P
            for k,e in mm:
                for _ in range(e):
                    t=G.pmul(t,expr[k])
                    if t==0: break
                if t==0: break
            if t!=0: out=G.padd(out,t)
        if not (isinstance(out,int) and out%P==0): res.append((a,out))
    return dict(nzc=nzc,incchecks=incchecks,res=res,rank=len(piv),nfree=n-len(piv),
                piv=piv,R=R,lin=lin,non=non,skipped=skipped)
if __name__=='__main__':
    for arg in sys.argv[1:]:
        FL=[int(x) for x in arg.split(',') if x] if arg!='-' else []
        v=list(v0)
        for b in FL: v[b]=1-v[b]
        ad.fwd(v,rounds=8)
        sc=L.NEQ-len(L.failing_eqs(L.all_atom_values(v)))
        r=analyse(v)
        tag='FULL MOD-P SOLUTION' if (not r['incchecks'] and not r['nzc'] and not r['res']) else ''
        print('flip %-22s score=%d rank=%d nfree=%d  ninc=%d %s | nzc=%d %s | res=%s  %s'
              %(arg,sc,r['rank'],r['nfree'],len(r['incchecks']),
                [(a,len(L.atom2eq.get(a,{}))) for a in r['incchecks']][:6],
                len(r['nzc']),[(a,len(L.atom2eq.get(a,{}))) for a,_ in r['nzc']][:6],
                [(a,('C' if isinstance(g,int) else 'P%d'%len(g))) for a,g in r['res']][:8], tag),flush=True)
