#!/usr/bin/env python3
"""Agent AB round 4: finish the solving-degree measurement.  Same tight (saturated) model as
ab_dreg2.py; two optimisations so n=4 fits in this box:
  (a) COLUMN SUPPORT restriction -- index only monomials that actually occur in some row (plus the
      target monomials).  All-zero columns cannot affect rank or row-space membership, so this is
      exact, not an approximation.
  (b) ONE augmented rank instead of n.  All n targets lie in the row space  <=>
      rank(M with the n target rows appended) == rank(M).  That is the only question asked.
"""
import sys, time, gc
from flint import nmod_mat
import ab_dreg2 as D

def solving_degree(n, dmax=7, cap_cells=4.0e8):
    dat=D.make(n)
    eqs,NV,BOOL,q,S=D.build_system(n,dat)
    print("  n=%d : %d vars (%d boolean), %d generators"%(n,NV,len(BOOL),len(eqs)),flush=True)
    for d in range(2,dmax+1):
        t0=time.time()
        rows=[]
        for p in eqs:
            dp=max(D.deg(e) for e in p)
            for m in D.monos_upto(NV,BOOL,d-dp):
                r={}; ok=True
                for e,c in p.items():
                    t=list(e)
                    for i in range(NV):
                        t[i]+=m[i]
                        if i in BOOL and t[i]>1: t[i]=1
                    tt=tuple(t)
                    if D.deg(tt)>d: ok=False; break
                    r[tt]=(r.get(tt,0)+c)%q
                if ok:
                    r={k:v for k,v in r.items() if v}
                    if r: rows.append(r)
        sup=set()
        for r in rows: sup.update(r)
        zero=tuple([0]*NV); sup.add(zero)
        tgts=[]
        for j in range(n):
            e=[0]*NV; e[j]=1; tt=tuple(e); sup.add(tt)
            tgts.append({tt:1, zero:(-S[j])%q})
        cols=sorted(sup); cidx={m:i for i,m in enumerate(cols)}
        cells=(len(rows)+n)*len(cols)
        print("     d=%d : %6d rows x %6d support-cols = %.2e cells (build %.0fs)"
              %(d,len(rows),len(cols),cells,time.time()-t0),flush=True)
        if cells>cap_cells:
            print("        OVER CAP (%.1e) -- cannot decide this degree on this box"%cap_cells,flush=True)
            return ('>=',d)
        M=nmod_mat(len(rows),len(cols),q)
        for i,r in enumerate(rows):
            for t,c in r.items(): M[i,cidx[t]]=c
        t1=time.time(); rk=M.rank(); del M; gc.collect()
        A=nmod_mat(len(rows)+n,len(cols),q)
        for i,r in enumerate(rows):
            for t,c in r.items(): A[i,cidx[t]]=c
        for j,tg in enumerate(tgts):
            for t,c in tg.items(): A[len(rows)+j,cidx[t]]=c
        rka=A.rank(); del A; gc.collect()
        print("        rank(M)=%d  rank(M+targets)=%d   %s   (%.0fs)"
              %(rk,rka,"ALL SELECTORS PINNED  <== SOLVING DEGREE" if rk==rka else "not yet",time.time()-t1),flush=True)
        if rk==rka: return ('=',d)
    return ('>',dmax)

if __name__=='__main__':
    lo,hi=int(sys.argv[1]),int(sys.argv[2])
    out={}
    for n in range(lo,hi+1):
        t0=time.time(); r=solving_degree(n); out[n]=r
        print("     -> n=%d solving degree %s%d   (%.0fs)\n"%(n,r[0],r[1],time.time()-t0),flush=True)
    print("MEASURED SOLVING DEGREES:", {k:("%s%d"%v) for k,v in out.items()})
