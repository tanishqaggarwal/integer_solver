#!/usr/bin/env python3
"""Agent AB: MEASURE the degree of regularity of the boolean-selector ECDLP system.

Singular is NOT installed in this container (checked: no binary anywhere on the filesystem),
so this measures d_reg directly with sympy's Groebner engine over GF(q).

CRITICAL MODELLING POINT (I nearly got this wrong):
the sibling MUST be the CURVE/ladder system, NOT the boolean modular subset-sum
  Sum s_i 2^i = k0 (mod q).
The subset-sum version is trivial once k0 is known -- that is §3's whole finding -- so it would
report a misleadingly low d_reg.  The unknown is k0 itself, so the ladder over the curve is the
only faithful sibling.

System, per step j (P_j = (a_j,b_j) a KNOWN point, s_j the selector, lam_j the chord slope):
  E1 : s*(lam*(x-a) - (y-b))                 = 0    (chord holds when the step is taken)
  E2 : x' - x - s*(lam^2 - 2x - a)           = 0
  E3 : y' - y - s*(lam*(x-x') - 2y)          = 0
  E4 : s^2 - s                               = 0
  E5 : (1-s)*lam                             = 0    (pins lam=0 on the untaken branch,
                                                     killing the spurious component)
endpoints (x_0,y_0)=R0 and (x_n,y_n)=R0+kG are CONSTANTS, so k is determined by the s_j only.
"""
import sys, time, random
from sympy import symbols, groebner, GF, Poly

def curve_setup(q,A,B):
    def add(P,Q):
        if P is None: return Q
        if Q is None: return P
        x1,y1=P; x2,y2=Q
        if (x1-x2)%q==0:
            if (y1+y2)%q==0: return None
            l=(3*x1*x1+A)*pow(2*y1%q,q-2,q)%q
        else:
            l=(y2-y1)*pow((x2-x1)%q,q-2,q)%q
        x3=(l*l-x1-x2)%q
        return (x3,(l*(x1-x3)-y1)%q)
    return add

def build(n, q=10007, A=2, B=3, seed=1):
    add=curve_setup(q,A,B)
    rng=random.Random(seed)
    pts=[(x,y) for x in range(q) for y in range(50) if (y*y-x*x*x-A*x-B)%q==0]
    G=None
    for x in range(1,q):
        r=(x*x*x+A*x+B)%q
        y=pow(r,(q+1)//4,q) if q%4==3 else None
        if y is not None and (y*y-r)%q==0: G=(x,y); break
    assert G is not None
    P=[G]
    for _ in range(n-1): P.append(add(P[-1],P[-1]))
    if any(p is None for p in P): return None
    for tries in range(400):
        k=rng.randrange(1,2**n)
        S=[(k>>j)&1 for j in range(n)]
        R=[G]                                # R0 = G  (offset, keeps us off infinity)
        ok=True
        for j in range(n):
            if S[j]:
                if R[-1][0]==P[j][0]: ok=False; break   # doubling / infinity: exceptional
                R.append(add(R[-1],P[j]))
            else:
                R.append(R[-1])
            if R[-1] is None: ok=False; break
        if ok and len(R)==n+1: return q,A,B,P,R,S,k
    return None

def system(n, dat):
    q,A,B,P,R,S,k=dat
    eqs=[]; names=[]
    sv=symbols('s0:%d'%n); lv=symbols('l0:%d'%n)
    xv=[None]*(n+1); yv=[None]*(n+1)
    xv[0],yv[0]=R[0]; xv[n],yv[n]=R[n]
    xi=symbols('x1:%d'%n); yi=symbols('y1:%d'%n)
    for j in range(1,n): xv[j],yv[j]=xi[j-1],yi[j-1]
    for j in range(n):
        a,b=P[j]; s=sv[j]; l=lv[j]; x=xv[j]; y=yv[j]; x2=xv[j+1]; y2=yv[j+1]
        eqs.append(s*(l*(x-a)-(y-b)))
        eqs.append(x2-x-s*(l*l-2*x-a))
        eqs.append(y2-y-s*(l*(x-x2)-2*y))
        eqs.append(s*s-s)
        eqs.append((1-s)*l)
    gens=list(sv)+list(lv)+list(xi)+list(yi)
    return eqs,gens,S,q

if __name__=='__main__':
    print("Singular present:", False, " -> measuring with sympy's Groebner engine over GF(q)")
    print()
    print("  n   vars  eqns   time     max deg of reduced degrevlex GB    #GB elts   recovers s?")
    for n in range(2, int(sys.argv[1]) if len(sys.argv)>1 else 6):
        dat=build(n)
        if dat is None: print("  %2d   (no clean instance)"%n); continue
        eqs,gens,S,q=system(n,dat)
        t0=time.time()
        try:
            G=groebner(eqs, *gens, order='grevlex', modulus=q)
        except Exception as e:
            print("  %2d   FAILED: %s"%(n,e)); continue
        el=time.time()-t0
        polys=[Poly(g,*gens,modulus=q) for g in G.exprs]
        md=max(p.total_degree() for p in polys)
        # does the GB pin every selector to its true value?
        pin=0
        for j,s in enumerate(S):
            tgt=gens[j]-s
            if any(Poly(tgt,*gens,modulus=q)==p or Poly(-tgt,*gens,modulus=q)==p for p in polys): pin+=1
        print("  %2d   %4d  %4d  %6.1fs            %2d                       %4d       %d/%d"
              %(n,len(gens),len(eqs),el,md,len(polys),pin,n))
        sys.stdout.flush()
