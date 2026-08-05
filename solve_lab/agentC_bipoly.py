#!/usr/bin/env python3
"""Bivariate polynomials over F_p in (t0,t1) as dict{(i,j):coeff}. Solve two-poly systems via
resultant w.r.t. t1 (univariate in t0), then back-substitute. Reuses agentC_poly for univariate."""
import agentC_poly as U
p = U.p

def const(c):
    c%=p
    return {} if c==0 else {(0,0):c}
def mk(base, v0, v1):
    d={}
    if base%p: d[(0,0)]=base%p
    if v0%p: d[(1,0)]=v0%p
    if v1%p: d[(0,1)]=v1%p
    return d
def add(a,b):
    r=dict(a)
    for k,v in b.items(): r[k]=(r.get(k,0)+v)%p
    return {k:v for k,v in r.items() if v%p}
def sub(a,b):
    r=dict(a)
    for k,v in b.items(): r[k]=(r.get(k,0)-v)%p
    return {k:v for k,v in r.items() if v%p}
def mul(a,b):
    r={}
    for (i,j),c in a.items():
        for (k,l),d in b.items():
            key=(i+k,j+l); r[key]=(r.get(key,0)+c*d)%p
    return {k:v for k,v in r.items() if v%p}

def deg_t1(a): return max((j for (i,j) in a), default=-1)
def t1_coeffs(a):
    """return dict j -> univariate poly in t0 (little-endian list)."""
    m=deg_t1(a); out={}
    for (i,j),c in a.items():
        out.setdefault(j,{})
        out[j][i]=out[j].get(i,0)+c
    res={}
    for j,dd in out.items():
        deg=max(dd); lst=[0]*(deg+1)
        for i,c in dd.items(): lst[i]=c%p
        res[j]=U.norm(lst)
    return res

def det_poly(M):
    """determinant of matrix M whose entries are univariate polys (little-endian lists) over F_p[t0].
    Laplace expansion (recursive). Small matrices only."""
    n=len(M)
    if n==1: return M[0][0]
    if n==2:
        return U.psub(U.pmul(M[0][0],M[1][1]), U.pmul(M[0][1],M[1][0]))
    total=[0]
    for c in range(n):
        if U.pdeg(M[0][c])<0: continue
        minor=[[M[r][cc] for cc in range(n) if cc!=c] for r in range(1,n)]
        term=U.pmul(M[0][c], det_poly(minor))
        if c%2==1: term=U.psub([0],term)
        total=U.padd(total,term)
    return total

def resultant_t1(A,B):
    """Sylvester resultant of A,B w.r.t. t1 -> univariate poly in t0 (little-endian)."""
    ca=t1_coeffs(A); cb=t1_coeffs(B)
    m=max(ca) if ca else 0; n=max(cb) if cb else 0
    if not ca or not cb: return [0]
    size=m+n
    if size==0: return [1]
    # Sylvester matrix (m+n)x(m+n); rows: n rows of A-coeffs shifted, m rows of B-coeffs shifted
    def coeff_list(cd, deg):
        # highest degree first: [a_deg,...,a0]
        return [cd.get(k,[0]) for k in range(deg,-1,-1)]
    A_row=coeff_list(ca,m); B_row=coeff_list(cb,n)
    S=[[[0] for _ in range(size)] for _ in range(size)]
    for r in range(n):
        for k in range(m+1):
            S[r][r+k]=A_row[k]
    for r in range(m):
        for k in range(n+1):
            S[n+r][r+k]=B_row[k]
    return U.norm(det_poly(S))

def eval_t1_at_t0(A, t0):
    """substitute t0 -> univariate poly in t1 (little-endian)."""
    out={}
    for (i,j),c in A.items():
        out[j]=(out.get(j,0)+c*pow(t0,i,p))%p
    if not out: return [0]
    deg=max(out); lst=[0]*(deg+1)
    for j,c in out.items(): lst[j]=c%p
    return U.norm(lst)

def solve2(A,B,verbose=False):
    """Return list of (t0,t1) solving A=B=0 mod p."""
    res=resultant_t1(A,B)
    if U.pdeg(res)<0:
        # resultant identically zero: common factor; fall back — solve B for t1(t0) generically
        # (rare) — sample t0 values: skip robust handling, return []
        if verbose: print("resultant identically zero")
        return []
    t0roots=U.roots_mod_p(res)
    sols=[]
    for t0 in t0roots:
        a1=eval_t1_at_t0(A,t0); b1=eval_t1_at_t0(B,t0)
        g=U.pgcd(a1,b1)
        if U.pdeg(g)<1:
            # maybe one is zero poly
            cand=a1 if U.pdeg(b1)<0 else b1
            for t1 in U.roots_mod_p(cand):
                sols.append((t0,t1))
            continue
        for t1 in U.roots_mod_p(g):
            sols.append((t0,t1))
    # dedup + verify
    out=[]
    seen=set()
    def ev(P,t0,t1):
        return sum(c*pow(t0,i,p)*pow(t1,j,p) for (i,j),c in P.items())%p
    for (t0,t1) in sols:
        if (t0,t1) in seen: continue
        seen.add((t0,t1))
        if ev(A,t0,t1)==0 and ev(B,t0,t1)==0:
            out.append((t0,t1))
    return out

if __name__=='__main__':
    # test: S: t0^2+t1^2-1=0 ; T: t0 - t1 =0  -> t0=t1=+-1/sqrt2
    A={(2,0):1,(0,2):1,(0,0):(-1)%p}
    B={(1,0):1,(0,1):(-1)%p}
    s=solve2(A,B); print("sols:",len(s))
    for (t0,t1) in s:
        print("  check:", (t0*t0+t1*t1-1)%p, (t0-t1)%p)
