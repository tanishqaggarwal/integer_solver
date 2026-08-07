"""Compute the closed non-boolean symbol set for a given state."""
import os, sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
import suppfree
_BC={}
def isbool(u):
    if u in _BC: return _BC[u]
    r=False
    for a in L.var_atoms[u]:
        pl=L.polys[a]
        if len(pl)==2 and (u,) in pl and (u,u) in pl and pl[(u,)]==-pl[(u,u)]: r=True;break
    _BC[u]=r; return r
def closure(v, seed=None, maxit=15, verbose=False):
    vm=[x%P for x in v]
    idx, freelist, vs = suppfree.build(vm, modp=True)
    FREESET=set(freelist)
    def supp(a):
        m=0
        for w in L.avars[a]: m |= vs[w] if w<len(vs) else 0
        s={freelist[i] for i in range(len(freelist)) if (m>>i)&1}
        s |= {w for w in L.avars[a] if w in FREESET}
        return s
    S=set(seed) if seed else set()
    if not S:
        # seed: non-boolean free inputs of every currently-nonzero check
        av=L.all_atom_values(v)
        for a in gsym.check_atoms():
            if av[a]%P:
                S |= {u for u in supp(a) if not isbool(u)}
    for it in range(maxit):
        SY=sorted(S); n=len(SY)
        val=gsym.build(v,SY,cap=None,verbose=False)
        sym=[]
        for a in gsym.check_atoms():
            f=gsym.evalpoly_sym(a,val,n,None)
            if not isinstance(f,int) or f%P: sym.append(a)
        new=set()
        for a in sym:
            new |= {u for u in supp(a) if not isbool(u)}
        if verbose: print('   closure it%d |S|=%d sym=%d new=%d'%(it,len(S),len(sym),len(new-S)))
        if new<=S: break
        S|=new
    return sorted(S)
