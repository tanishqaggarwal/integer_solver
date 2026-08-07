#!/usr/bin/env python3
"""Agent AB: measure d_reg properly.  Two fixes over ab_dreg.py:

FIX 1 -- the model was LEAKY.  With only E1..E5, whenever an intermediate point R_j coincides
with the point P_j being added, E1 reads s*(lam*0 - 0) = 0 and lam is UNCONSTRAINED: a spurious
POSITIVE-DIMENSIONAL component.  That is why the n>=3 runs returned max-GB-degree 2 and failed to
pin the selectors.  Saturate it away (Rabinowitsch):
     E6 : u*(x-a) - s = 0      (s=1  =>  x != a and u = 1/(x-a))
     E7 : (1-s)*u   = 0        (s=0  =>  u = 0)

FIX 2 -- max degree of the REDUCED GB is the wrong statistic.  For a zero-dimensional ideal with a
unique solution the reduced GB is {v - c_v}, degree 1, for every n.  The quantity in my complexity
claim is the SOLVING DEGREE: the smallest d at which the degree-d Macaulay row space already
contains (v - c_v).  Measured here by explicit Macaulay construction + rref over GF(q) (flint),
working in the quotient by s_j^2 = s_j so the s-part of every monomial is multilinear.
"""
import sys, time, random, itertools
from flint import nmod_mat

def make(n, q=10007, A=2, B=3, seed=1):
    def add(P,Q):
        if P is None: return Q
        if Q is None: return P
        x1,y1=P; x2,y2=Q
        if (x1-x2)%q==0:
            if (y1+y2)%q==0: return None
            l=(3*x1*x1+A)*pow(2*y1%q,q-2,q)%q
        else: l=(y2-y1)*pow((x2-x1)%q,q-2,q)%q
        x3=(l*l-x1-x2)%q
        return (x3,(l*(x1-x3)-y1)%q)
    G=None
    for x in range(1,q):
        r=(x*x*x+A*x+B)%q
        y=pow(r,(q+1)//4,q)
        if (y*y-r)%q==0: G=(x,y); break
    P=[G]
    for _ in range(n-1): P.append(add(P[-1],P[-1]))
    rng=random.Random(seed)
    for _ in range(500):
        k=rng.randrange(1,2**n); S=[(k>>j)&1 for j in range(n)]
        R=[G]; ok=True
        for j in range(n):
            if S[j]:
                if R[-1][0]==P[j][0]: ok=False;break
                R.append(add(R[-1],P[j]))
                if R[-1] is None: ok=False;break
            else: R.append(R[-1])
        if ok: return q,P,R,S
    return None

# ---- polynomials as dict{ exponent-tuple : coeff }, quotient by s_j^2=s_j -----------------
def build_system(n, dat):
    q,P,R,S=dat
    # variable order: s_0..s_{n-1}, l_0.., u_0.., x_1..x_{n-1}, y_1..y_{n-1}
    NV=3*n+2*(n-1)
    def si(j): return j
    def li(j): return n+j
    def ui(j): return 2*n+j
    def xi(j): return 3*n+(j-1)          # j=1..n-1
    def yi(j): return 3*n+(n-1)+(j-1)
    BOOL=set(range(n))
    def mono(**kw):
        e=[0]*NV
        for k_,v in kw.items(): pass
        return tuple(e)
    def V(i):
        e=[0]*NV; e[i]=1; return {tuple(e):1}
    def C(c): return {tuple([0]*NV): c%q} if c%q else {}
    def mul(a,b):
        out={}
        for ea,ca in a.items():
            for eb,cb in b.items():
                e=list(ea)
                ok=True
                for i in range(NV):
                    e[i]+=eb[i]
                    if i in BOOL and e[i]>1: e[i]=1     # s^2 = s
                t=tuple(e); v=(out.get(t,0)+ca*cb)%q
                if v: out[t]=v
                elif t in out: del out[t]
        return out
    def addp(*ps):
        out={}
        for p in ps:
            for e,c in p.items():
                v=(out.get(e,0)+c)%q
                if v: out[e]=v
                elif e in out: del out[e]
        return out
    def neg(p): return {e:(-c)%q for e,c in p.items()}
    eqs=[]
    for j in range(n):
        a,b=P[j]
        s=V(si(j)); l=V(li(j)); u=V(ui(j))
        x = V(xi(j)) if 1<=j<=n-1 else C(R[j][0])
        y = V(yi(j)) if 1<=j<=n-1 else C(R[j][1])
        x2= V(xi(j+1)) if 1<=j+1<=n-1 else C(R[j+1][0])
        y2= V(yi(j+1)) if 1<=j+1<=n-1 else C(R[j+1][1])
        one=C(1)
        eqs.append(mul(s, addp(mul(l, addp(x,C(-a))), neg(addp(y,C(-b))))))          # E1
        eqs.append(addp(x2, neg(x), neg(mul(s, addp(mul(l,l), neg(mul(C(2),x)), C(-a))))))  # E2
        eqs.append(addp(y2, neg(y), neg(mul(s, addp(mul(l,addp(x,neg(x2))), neg(mul(C(2),y)))))))  # E3
        eqs.append(addp(mul(s,s), neg(s)))                                            # E4 (=0 in quotient)
        eqs.append(mul(addp(one,neg(s)), l))                                          # E5
        eqs.append(addp(mul(u, addp(x,C(-a))), neg(s)))                               # E6 saturation
        eqs.append(mul(addp(one,neg(s)), u))                                          # E7
    eqs=[e for e in eqs if e]
    return eqs, NV, BOOL, q, S

def monos_upto(NV,BOOL,d):
    """all exponent tuples of total degree <= d, s-part multilinear"""
    out=[]
    def rec(i,rem,cur):
        if i==NV:
            out.append(tuple(cur)); return
        hi = 1 if i in BOOL else rem
        for e in range(0,min(hi,rem)+1):
            cur.append(e); rec(i+1,rem-e,cur); cur.pop()
    rec(0,d,[])
    return out

def deg(e): return sum(e)

def solving_degree(n, dmax=7, cap_cols=200000, cap_cells=6*10**8):
    dat=make(n)
    eqs,NV,BOOL,q,S=build_system(n,dat)
    print("  n=%d : %d vars (%d boolean), %d generators, max gen degree %d"
          %(n,NV,len(BOOL),len(eqs),max(max(deg(e) for e in p) for p in eqs)))
    for d in range(2,dmax+1):
        cols=monos_upto(NV,BOOL,d)
        cidx={m:i for i,m in enumerate(cols)}
        rows=[]
        for p in eqs:
            dp=max(deg(e) for e in p)
            for m in monos_upto(NV,BOOL,d-dp):
                r={}
                for e,c in p.items():
                    t=list(e); ok=True
                    for i in range(NV):
                        t[i]+=m[i]
                        if i in BOOL and t[i]>1: t[i]=1
                    tt=tuple(t)
                    if deg(tt)>d: ok=False;break
                    r[tt]=(r.get(tt,0)+c)%q
                if ok and any(r.values()): rows.append(r)
        cells=len(rows)*len(cols)
        if len(cols)>cap_cols or cells>cap_cells:
            print("     d=%d : %d rows x %d cols = %.2e cells -- OVER CAP, stopping"%(d,len(rows),len(cols),cells))
            return None
        M=nmod_mat(len(rows),len(cols),q)
        for i,r in enumerate(rows):
            for t,c in r.items():
                if c: M[i,cidx[t]]=c
        t0=time.time(); R=M.rref()[0]; el=time.time()-t0
        rk=M.rank()
        # extract pivots of the rref, then reduce each target vector against them (cheap)
        piv=[]; used=set()
        for i in range(rk):
            for jj in range(len(cols)):
                if R[i,jj]!=0:
                    piv.append((i,jj)); used.add(jj); break
        def in_rowspace(vec):
            v=dict(vec)
            for i,jj in piv:
                c=v.get(jj,0)
                if c:
                    inv=pow(int(R[i,jj]),q-2,q)
                    f=c*inv%q
                    for kk in range(len(cols)):
                        rv=int(R[i,kk])
                        if rv:
                            nv=(v.get(kk,0)-f*rv)%q
                            if nv: v[kk]=nv
                            elif kk in v: del v[kk]
            return not any(v.values())
        pinned=0
        for j in range(n):
            e=[0]*NV; e[j]=1
            vec={cidx[tuple(e)]:1, cidx[tuple([0]*NV)]:(-S[j])%q}
            if in_rowspace(vec): pinned+=1
        print("     d=%d : %6d rows x %6d cols  rank %6d  rref %5.1fs  selectors pinned %d/%d %s"
              %(d,len(rows),len(cols),rk,el,pinned,n,"  <== SOLVING DEGREE" if pinned==n else ""))
        sys.stdout.flush()
        if pinned==n: return d
    return None

if __name__=='__main__':
    print("MEASURED SOLVING DEGREE of the boolean-selector ECDLP ladder (tight/saturated model)")
    print()
    res={}
    for n in range(2, int(sys.argv[1]) if len(sys.argv)>1 else 5):
        t0=time.time(); d=solving_degree(n); res[n]=d
        print("     -> n=%d solving degree = %s   (%.1fs)"%(n,d,time.time()-t0)); print()
    print("MEASURED:", res)
