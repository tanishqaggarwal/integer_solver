"""Exact symbolic forward evaluation over F_p, SPARSE monomials.
monomial = tuple of (symbol_index, exponent) sorted by index.  Scales to thousands
of symbols.  Polynomial = int (constant) or dict{monomial: coeff}.
"""
import os, sys, time
LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff'))
sys.path.insert(0, os.path.join(LAB,'s10'))
import lib as L, ad
P = ad.P
EMPTY = ()

def mmul(m1, m2, boolset=None):
    if not m1: return m2
    if not m2: return m1
    d={}
    for k,e in m1: d[k]=e
    for k,e in m2: d[k]=d.get(k,0)+e
    if boolset:
        for k in list(d):
            if k in boolset and d[k]>1: d[k]=1
    return tuple(sorted(d.items()))

def mdeg(m): return sum(e for _,e in m)

def pmul(f,g,cap=None,boolset=None):
    if isinstance(f,int) and isinstance(g,int): return f*g%P
    if isinstance(f,int):
        if f==0: return 0
        return {m:(c*f)%P for m,c in g.items() if (c*f)%P}
    if isinstance(g,int):
        if g==0: return 0
        return {m:(c*g)%P for m,c in f.items() if (c*g)%P}
    out={}
    for m1,c1 in f.items():
        for m2,c2 in g.items():
            m=mmul(m1,m2,boolset)
            if cap and mdeg(m)>cap: raise OverflowError('cap')
            out[m]=(out.get(m,0)+c1*c2)%P
    out={m:c for m,c in out.items() if c}
    if not out: return 0
    if len(out)==1 and EMPTY in out: return out[EMPTY]
    return out

def padd(f,g):
    if isinstance(f,int) and isinstance(g,int): return (f+g)%P
    if isinstance(f,int):
        if f%P==0: return g
        f={EMPTY:f%P}
    if isinstance(g,int):
        if g%P==0: return f
        g={EMPTY:g%P}
    out=dict(f)
    for m,c in g.items():
        c2=(out.get(m,0)+c)%P
        if c2: out[m]=c2
        elif m in out: del out[m]
    if not out: return 0
    if len(out)==1 and EMPTY in out: return out[EMPTY]
    return out

def deg(f): return 0 if isinstance(f,int) else max(mdeg(m) for m in f)
def nterms(f): return 1 if isinstance(f,int) else len(f)

def evalatom(a, val, cap=None, skip=None, boolset=None):
    s=0
    for m,c in L.polys[a].items():
        if skip is not None and skip in m: continue
        t=c%P
        for u in m:
            t=pmul(t,val[u],cap,boolset)
            if t==0: break
        if t!=0: s=padd(s,t)
    return s

def build(v, syms, cap=None, boolsyms=(), verbose=False):
    sidx={u:i for i,u in enumerate(syms)}
    boolset={sidx[u] for u in boolsyms if u in sidx}
    val=[x%P for x in v]
    for u in syms: val[u]={((sidx[u],1),):1}
    symset=set(syms)
    skipped=0
    t0=time.time()
    for t in ad.ORDER:
        if t in symset: continue
        a=L.definer[t]
        c=0; ok=True
        for m,cc in L.polys[a].items():
            k=m.count(t)
            if k>1: ok=False; break
            if k==1:
                tt=cc%P
                for u in m:
                    if u!=t: tt=pmul(tt,val[u],cap,boolset)
                c=padd(c,tt)
        if not ok or c==0 or not isinstance(c,int):
            if not isinstance(c,int): skipped+=1
            continue
        rest=evalatom(a,val,cap,skip=t,boolset=boolset)
        val[t]=pmul(rest,(-pow(c,-1,P))%P,cap,boolset)
    if verbose:
        print('   sym pass %.1fs skipped_nonlinear_gates=%d'%(time.time()-t0,skipped))
    return val, skipped

def check_atoms():
    return [a for a in range(L.NA) if a not in L.atom_out]
