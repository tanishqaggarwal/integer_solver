"""Exact equation-level mod-p ceiling of a FRAME = (boolean flips, detached gate vars)."""
import os, sys, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G, gsolve, gGclose
from gsym2 import L, ad, P

def build_detached(v, DETACH, syms, cap=6):
    """Symbolic forward eval with DETACH variables treated as free symbols."""
    DET=set(DETACH)
    sidx={u:i for i,u in enumerate(syms)}
    val=[x%P for x in v]
    for u in syms: val[u]={((sidx[u],1),):1}
    symset=set(syms)
    for t in ad.ORDER:
        if t in symset or t in DET: continue
        a=L.definer[t]
        c=0; ok=True
        for m,cc in L.polys[a].items():
            k=m.count(t)
            if k>1: ok=False; break
            if k==1:
                tt=cc%P
                for u in m:
                    if u!=t: tt=G.pmul(tt,val[u],cap)
                c=G.padd(c,tt)
        if not ok or c==0 or not isinstance(c,int): continue
        rest=G.evalatom(a,val,cap,skip=t)
        val[t]=G.pmul(rest,(-pow(c,-1,P))%P,cap)
    return val

def ceiling(v, DETACH=(), cap=6, verbose=False, want_res=False):
    DET=set(DETACH)
    FREE=[u for u in range(L.NVARS) if u not in L.definer or u in DET]
    NB=[u for u in FREE if not gGclose.isbool(u)]
    n=len(NB)
    val=build_detached(v,DET,NB,cap)
    allsym=[G.evalatom(a,val,cap) for a in range(L.NA)]
    bad=[]
    for i,(m,sq,co) in enumerate(L.eq_atoms):
        s=0
        for a,c in co.items():
            t=allsym[a]
            if isinstance(t,int):
                if t%P: s=G.padd(s,(c*t)%P)
            else: s=G.padd(s,{mm:(c*cc)%P for mm,cc in t.items()})
        if not (isinstance(s,int) and s%P==0): bad.append((i,s))
    lin=[(i,s) for i,s in bad if not isinstance(s,int) and G.deg(s)==1]
    non=[(i,s) for i,s in bad if isinstance(s,int) or G.deg(s)>1]
    sp=[]
    for i,s in lin:
        r={}
        for mm,c in s.items():
            if not mm: r[n]=(-c)%P
            else: r[mm[0][0]]=c%P
        sp.append(r)
    piv,R=gsolve.sparse_rref(sp,n)
    inc=[lin[i][0] for i,r in enumerate(R) if not [c for c in r if c!=n] and r.get(n,0)%P]
    need=set()
    for i,s in non:
        if isinstance(s,int): continue
        for mm in s:
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
    for i,s in non:
        if isinstance(s,int):
            if s%P: res.append((i,s%P))
            continue
        out=0
        for mm,c in s.items():
            t=c%P
            for k,e in mm:
                for _ in range(e):
                    t=G.pmul(t,expr[k])
                    if t==0: break
                if t==0: break
            if t!=0: out=G.padd(out,t)
        if not (isinstance(out,int) and out%P==0): res.append((i,out))
    nfail=len(inc)+len(res)
    if want_res: return nfail,inc,res,NB,piv,R,n
    return nfail,inc,res
