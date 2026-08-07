"""Exact symbolic forward evaluation of the circuit over F_p.

Gate output coefficients are +-1, so forward evaluation is an honest polynomial map
from the free inputs to every variable.  Choosing a set S of free inputs as
indeterminates and evaluating everything else at a consistent state gives, for every
CHECK atom, an exact polynomial in F_p[S].  That is the reduced system.
"""
import os, sys, time, collections
LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff'))
sys.path.insert(0, os.path.join(LAB,'s10'))
import lib as L, ad
P = ad.P

# ---------- sparse multivariate polys over F_p, monomial = tuple of exponents ----
class Ctx:
    def __init__(self, n): self.n = n; self.ZERO=(0,)*n

def pmul(f, g, n, cap):
    if isinstance(f,int) and isinstance(g,int): return f*g % P
    if isinstance(f,int):
        if f==0: return 0
        return {m:(c*f)%P for m,c in g.items() if (c*f)%P}
    if isinstance(g,int):
        if g==0: return 0
        return {m:(c*g)%P for m,c in f.items() if (c*g)%P}
    out={}
    for m1,c1 in f.items():
        for m2,c2 in g.items():
            m=tuple(a+b for a,b in zip(m1,m2))
            if cap and sum(m)>cap: raise OverflowError('degree cap')
            c=out.get(m,0)+c1*c2
            out[m]=c%P
    out={m:c for m,c in out.items() if c}
    if not out: return 0
    if len(out)==1 and next(iter(out))==(0,)*n: return out[(0,)*n]
    return out

def padd(f,g,n):
    if isinstance(f,int) and isinstance(g,int): return (f+g)%P
    if isinstance(f,int):
        if f==0: return g
        f={(0,)*n: f%P} if f%P else {}
    if isinstance(g,int):
        if g==0: return f
        g={(0,)*n: g%P} if g%P else {}
    out=dict(f)
    for m,c in g.items():
        c2=(out.get(m,0)+c)%P
        if c2: out[m]=c2
        elif m in out: del out[m]
    if not out: return 0
    if len(out)==1 and next(iter(out))==(0,)*n: return out[(0,)*n]
    return out

def deg(f):
    if isinstance(f,int): return 0
    return max(sum(m) for m in f)

def nterms(f): return 1 if isinstance(f,int) else len(f)

def evalpoly_sym(a, val, n, cap, skip=None):
    """symbolic value of atom a; val[u] is int-or-dict; skip: variable to omit."""
    s = 0
    for m,c in L.polys[a].items():
        if skip is not None and skip in m: continue
        t = c % P
        for u in m:
            t = pmul(t, val[u], n, cap)
            if t == 0: break
        if t != 0:
            s = padd(s, t, n)
    return s

def build(v, syms, cap=12, verbose=True):
    """v: integer state (already forward-consistent). syms: list of free input indices.
    Returns val[] list of symbolic values (int or dict) for every variable."""
    n = len(syms)
    sidx = {u:i for i,u in enumerate(syms)}
    val = [x % P for x in v]
    for u in syms:
        e=[0]*n; e[sidx[u]]=1
        val[u] = {tuple(e):1}
    symset=set(syms)
    t0=time.time(); maxd=0; maxt=0
    for t in ad.ORDER:
        if t in symset: continue
        a = L.definer[t]
        # coefficient of t in atom a (must be linear in t)
        c = 0; ok=True
        for m,cc in L.polys[a].items():
            k=m.count(t)
            if k>1: ok=False; break
            if k==1:
                # coefficient may itself be symbolic -> require it constant
                tt=cc%P
                for u in m:
                    if u!=t: tt=pmul(tt,val[u],n,cap)
                c = padd(c,tt,n)
        if not ok or c==0:
            continue
        if not isinstance(c,int):
            # nonlinear gate: skip (leave numeric value)
            continue
        rest = evalpoly_sym(a, val, n, cap, skip=t)
        inv = pow(c,-1,P)
        val[t] = pmul(rest, (-inv)%P, n, cap)
        d=deg(val[t]); maxd=max(maxd,d); maxt=max(maxt,nterms(val[t]))
    if verbose:
        print('  symbolic pass: %.1fs  maxdeg=%d maxterms=%d' % (time.time()-t0,maxd,maxt))
    return val

def check_atoms():
    return [a for a in range(L.NA) if a not in L.atom_out]
