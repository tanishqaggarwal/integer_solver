"""Sparse linear solve over F_p + substitution into the nonlinear part."""
import sys, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import P

def sparse_rref(rows, ncol, verbose=False):
    """rows: list of dict {col:coef} with col==ncol meaning RHS(const moved to rhs).
    Returns (pivots dict col->row, rows) in reduced form."""
    R=[dict(r) for r in rows]
    piv={}
    order=sorted(range(len(R)), key=lambda i: len(R[i]))
    t0=time.time()
    for it,i in enumerate(order):
        r=R[i]
        # reduce by existing pivots
        changed=True
        while changed:
            changed=False
            for c in list(r.keys()):
                if c!=ncol and c in piv:
                    f=r[c]; pr=R[piv[c]]
                    for c2,v2 in pr.items():
                        nv=(r.get(c2,0)-f*v2)%P
                        if nv: r[c2]=nv
                        elif c2 in r: del r[c2]
                    changed=True
                    break
        cols=[c for c in r if c!=ncol]
        if not cols:
            continue
        c=min(cols, key=lambda c:(len(r),c))
        iv=pow(r[c],-1,P)
        for k in list(r): r[k]=r[k]*iv%P
        piv[c]=i
        if verbose and it%200==0: print('    rref %d/%d piv=%d %.0fs'%(it,len(R),len(piv),time.time()-t0),flush=True)
    # back-substitute so every pivot row is reduced w.r.t. other pivots
    for c,i in list(piv.items()):
        r=R[i]
        changed=True
        while changed:
            changed=False
            for c2 in list(r.keys()):
                if c2!=ncol and c2!=c and c2 in piv:
                    f=r[c2]; pr=R[piv[c2]]
                    for c3,v3 in pr.items():
                        nv=(r.get(c3,0)-f*v3)%P
                        if nv: r[c3]=nv
                        elif c3 in r: del r[c3]
                    changed=True; break
    return piv,R

def inconsistent(rows,ncol,piv,R):
    out=[]
    for i,r in enumerate(R):
        cols=[c for c in r if c!=ncol]
        if not cols and r.get(ncol,0)%P: out.append(i)
    return out
