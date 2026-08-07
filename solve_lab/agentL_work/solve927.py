"""Fit-and-solve for the c>1 divisibility conditions.  Method + guards inherited from P.
Never brute-force E.run over c; brute-force the FITTED polynomial instead (cheap), and always
verify by DIRECT RECOMPUTATION."""
import sys, pickle, collections, json, time
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentL_work/mkassign2.py').read().split('#MAINSTART')[0]
exec(src)
CGT={a for a,s in SL.items() if s and s%p==0 and abs(s)//p>1}
def factor(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def probe(vv,i,ws,ts):
    old=[vv[w] for w in ws]
    for w,t in zip(ws,ts): vv[w]+=p*t
    r=E.run(vv)[i]
    for w,o in zip(ws,old): vv[w]=o
    return r
def fit(vv,i,w,deg=4):
    """exact Newton forward-difference fit of R(t)/p as a polynomial in t; returns coeffs list"""
    ys=[]
    for t in range(deg+1):
        y=probe(vv,i,[w],[t])
        if y%p: return None,None
        ys.append(y//p)
    d=[ys[:]]
    for k in range(deg):
        d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
    coeffs=[d[k][0] for k in range(deg+1)]
    topdeg=max([k for k in range(deg+1) if coeffs[k]!=0], default=0)
    return coeffs, topdeg
def peval(coeffs,t,m):
    """Newton form: sum_k coeffs[k] * C(t,k), evaluated mod m"""
    tot=0; binom=1
    for k,c in enumerate(coeffs):
        if k>0: binom=binom*(t-k+1)//k
        tot=(tot+c*binom)%m
    return tot%m
def roots_mod(coeffs,q,e):
    m=q**e
    return [t for t in range(m) if peval(coeffs,t,m)==0]
def crt(pairs):
    t,mod=0,1
    for r,m in pairs:
        g=gcd(mod,m)
        if (r-t)%g: return None
        lcm=mod//g*m
        t=(t+mod*(((r-t)//g)*pow(mod//g,-1,m//g)%(m//g)))%lcm
        mod=lcm
    return t,mod
def solve_one(vv,a):
    """find a shift on some wire of atom a that clears its condition; exact, verified"""
    i=E.residx[a]; c=abs(SL[a])//p
    ws=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
    ws+=[q for q in atomvalvars[a] if q in SHIFT and q not in ws]
    F=factor(c)
    for w in ws:
        coeffs,td=fit(vv,i,w)
        if coeffs is None: continue
        if coeffs[0]%c==0: return (w,0)
        pairs=[]; ok=True
        for q,e in F.items():
            rs=roots_mod(coeffs,q,e)
            if not rs: ok=False; break
            pairs.append((rs[0],q**e))
        if not ok: continue
        r=crt(pairs)
        if r is None: continue
        t=r[0]
        # GUARD: verify by direct recomputation
        val=probe(vv,i,[w],[t])
        if val%(c*p)==0: return (w,t,td)
    return None
if __name__=='__main__':
    import random
    rnd=random.Random(7)
    for S in ([24601,2081], rnd.sample(M['live'],17)):
        v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
        vv=[0]*NV
        for k,x in v.items(): vv[k]=x
        for rd in range(60):   # greedy first, to its fixpoint
            bad=relift(vv)
            if not bad: break
            r=E.run(vv); fx=0
            for a in bad:
                i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
                if cur%p: continue
                imm=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
                for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                    old=vv[w]; vv[w]=old+p; d=E.run(vv)[i]-cur; vv[w]=old
                    if d==0: continue
                    g=gcd(d,sm)
                    if cur%g: continue
                    mm=sm//g
                    t=(-(cur//g))*pow((d//g)%mm,-1,mm)%mm if mm>1 else 0
                    vv[w]=old+p*t; fx+=1; break
            if fx==0: break
        print('\n=== |S|=%d : fit-and-solve on the greedy fixpoint ==='%len(S),flush=True)
        for outer in range(12):
            left=relift(vv); r=E.run(vv)
            stuck=[a for a in left if r[E.residx[a]]%p==0]
            if not stuck:
                print('  ALL c>1 CONDITIONS DISCHARGED after %d outer rounds'%outer,flush=True); break
            print('  round %d: %d stuck'%(outer,len(stuck)),flush=True)
            prog=0
            for a in stuck:
                t0=time.time(); res=solve_one(vv,a)
                if res and len(res)==3:
                    w,t,td=res; vv[w]+=p*t; prog+=1
                    print('     SOLVED c=%-9d deg=%d wire x%-6d t=%-10d (%.0fs, verified exactly)'%(
                        abs(SL[a])//p,td,w,t,time.time()-t0),flush=True)
                elif res: pass
                else: print('     no root: c=%d  %s'%(abs(SL[a])//p,a[:60]),flush=True)
            if prog==0: break
        left=relift(vv); r=E.run(vv)
        stuck=[a for a in left if r[E.residx[a]]%p==0]
        nz=sum(1 for x in r if x)
        print('  FINAL: %d undischarged, %d nonzero atoms (2 = the target congruences)'%(len(stuck),nz),flush=True)
