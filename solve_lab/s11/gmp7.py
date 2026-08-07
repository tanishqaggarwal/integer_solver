"""Solve the mod-p system by closure, with response-based knob discovery.

The cached matrix (gmp6) says which free input moves which check.  Grow an active set:
rows = the failing checks; knobs = every live input whose response touches a row; then add every
check those knobs disturb as a new row; repeat.  Solve over GF(p) with sparse elimination, apply,
re-evaluate, and Newton-iterate (responses for the active knobs are recomputed exactly at the new
base, which is cheap -- only the discovery matrix is expensive, and it is reused).

Booleans take care of themselves: perturbing a message bit breaks its own x^2-x check, that check
enters as a row, and the solver is then forced to give that knob coefficient 0.
"""
import sys, os, json, time, pickle, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]

def allchk(v): return {a:evalp(L.polys[a],v) for a in CHK}
def failing(d): return sorted(a for a,x in d.items() if x)

def resp_at(base, bd, u, rows=None):
    v=list(base); v[u]=(v[u]+1)%P; forwardp(v)
    d={}
    it = rows if rows is not None else CHK
    for a in it:
        x=evalp(L.polys[a],v)
        if x!=bd[a]: d[a]=(x-bd[a])%P
    return d

def sparse_solve(cols, rowl, rhs):
    """cols: list of {row_index: coeff}; solve sum_j x_j cols[j] = rhs over GF(p)."""
    n=len(cols); m=len(rowl)
    R=[dict() for _ in range(m)]
    for j,c in enumerate(cols):
        for i,x in c.items(): R[i][j]=x
    b=list(rhs)
    piv={}; used=[False]*m; order=[]
    for c in range(n):
        cand=[i for i in range(m) if not used[i] and R[i].get(c)]
        if not cand: continue
        i=min(cand, key=lambda i: len(R[i]))
        used[i]=True; piv[c]=i; order.append(c)
        inv=pow(R[i][c],-1,P)
        R[i]={k:x*inv%P for k,x in R[i].items()}; b[i]=b[i]*inv%P
        Ri=R[i]; bi=b[i]
        for k in range(m):
            if k==i: continue
            f=R[k].get(c)
            if not f: continue
            Rk=R[k]
            for kk,x in Ri.items():
                nv=(Rk.get(kk,0)-f*x)%P
                if nv: Rk[kk]=nv
                elif kk in Rk: del Rk[kk]
            b[k]=(b[k]-f*bi)%P
    for i in range(m):
        if not used[i] and b[i] and not R[i]: return None
    for i in range(m):
        if not R[i] and b[i]: return None
    x=[0]*n
    for c,i in piv.items(): x[c]=b[i]
    return x

def run(base, disc, maxnewton=8, maxgrow=6, verbose=True):
    for nt in range(maxnewton):
        bd=allchk(base); F=failing(bd)
        if verbose: print(f"NEWTON {nt}: failing checks mod p = {len(F)} {F[:16]}", flush=True)
        if not F: return base
        rows=set(F); knobs=[]
        for g in range(maxgrow):
            newk=[u for u,d in disc.items() if u not in knobs and (set(d) & rows)]
            if not newk and knobs: break
            knobs.extend(newk)
            R2=set(rows)
            for u in knobs: R2 |= set(disc.get(u,{}))
            if verbose: print(f"   grow {g}: knobs={len(knobs)} rows={len(R2)}", flush=True)
            if R2==rows and newk==[]: break
            rows=R2
            if len(rows)>2500 or len(knobs)>2500:
                if verbose: print("   closure too large, stopping growth"); break
        rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
        t0=time.time()
        cols=[]
        for u in knobs:
            d=resp_at(base,bd,u,rowl)
            cols.append({idx[a]:x for a,x in d.items()})
        if verbose: print(f"   exact responses at this base: {time.time()-t0:.0f}s", flush=True)
        rhs=[(-bd[a])%P for a in rowl]
        t0=time.time()
        x=sparse_solve(cols,rowl,rhs)
        if verbose:
            print(f"   system {len(rowl)}x{len(knobs)} -> {'SOLVED' if x is not None else 'INCONSISTENT'}"
                  f" ({time.time()-t0:.0f}s)", flush=True)
        if x is None: return None
        v=list(base)
        for j,u in enumerate(knobs): v[u]=(v[u]+x[j])%P
        forwardp(v)
        base=v
    return None

if __name__=='__main__':
    D=pickle.load(open(os.path.join(HERE,'data','resp_modp.pkl'),'rb'))
    base=D['base']; disc=D['cols']
    print(f"discovery matrix: {len(disc)} live knobs")
    t0=time.time()
    v=run(base, disc)
    print(f"total {time.time()-t0:.0f}s")
    if v is not None:
        json.dump([int(x) for x in v], open(os.path.join(HERE,'data','gmp7_state.json'),'w'))
        print("*** EVERY CHECK ZERO MOD P -- saved data/gmp7_state.json")
