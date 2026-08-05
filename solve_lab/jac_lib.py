"""Exact forward-tangent (dual number) first-order Jacobian over GF(p) at a loaded point.
Every gate is degree<=2; equations are polynomials. Dual numbers give exact d/du_j."""
import heal_harness as H
p=H.p
freelist=sorted(H.freeinp)
freeidx={j:k for k,j in enumerate(freelist)}   # free input var -> column index
NF=len(freelist)

class D:
    __slots__=('v','g')
    def __init__(self,v,g=None):
        self.v=v%p
        self.g=g if g is not None else {}
    def __add__(s,o):
        if isinstance(o,D):
            g=dict(s.g)
            for k,c in o.g.items():
                nc=(g.get(k,0)+c)%p
                if nc: g[k]=nc
                elif k in g: del g[k]
            return D((s.v+o.v)%p,g)
        return D((s.v+o)%p,dict(s.g))
    __radd__=__add__
    def __sub__(s,o):
        if isinstance(o,D):
            g=dict(s.g)
            for k,c in o.g.items():
                nc=(g.get(k,0)-c)%p
                if nc: g[k]=nc
                elif k in g: del g[k]
            return D((s.v-o.v)%p,g)
        return D((s.v-o)%p,dict(s.g))
    def __rsub__(s,o):  # o - s , o int
        g={k:(-c)%p for k,c in s.g.items()}
        return D((o-s.v)%p,g)
    def __mul__(s,o):
        if isinstance(o,D):
            g={}
            for k,c in s.g.items():
                nc=(c*o.v)%p
                if nc: g[k]=nc
            for k,c in o.g.items():
                nc=(g.get(k,0)+c*s.v)%p
                if nc: g[k]=nc
                elif k in g: del g[k]
            return D((s.v*o.v)%p,g)
        else:
            g={}
            for k,c in s.g.items():
                nc=(c*o)%p
                if nc: g[k]=nc
            return D((s.v*o)%p,g)
    __rmul__=__mul__
    def __neg__(s):
        return D((-s.v)%p,{k:(-c)%p for k,c in s.g.items()})

def build_duals():
    """Return list vd of DualSparse for every variable at current H.val."""
    val=H.val
    vd=[None]*H.NVARS
    for j in H.freeinp:
        vd[j]=D(val[j],{freeidx[j]:1})
    ns={'v':vd,'__builtins__':{}}
    for k,t in enumerate(H.order):
        r=eval(H.gcode[k],ns)
        vd[t]=r if isinstance(r,D) else D(r)
    return vd

def eq_jac_row(i,vd):
    """Return (residual_val mod p, grad dict) for equation i."""
    r=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
    if isinstance(r,D): return r.v, r.g
    return r%p, {}
