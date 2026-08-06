"""Mod-p closure solve.

Mod p a knob's response is affine and there is no divisibility test to fail, so the repair is a
LINEAR system over GF(p) rather than the integer cascade that stalled all session.  Build it by
closure: start from the failing checks, take the free inputs that move them, add every check
those inputs disturb as a new row, and repeat until the system closes.  Then solve.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp, solvep
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=[u for u in range(L.NVARS) if u not in L.definer]
FS=set(FREE)

def allchk(v): return {a:evalp(L.polys[a],v) for a in CHK}
def failing(d): return sorted(a for a,x in d.items() if x)

_cache={}
def response(base, bd, u):
    """exact affine response of every check to u -> u+1, as {check: delta}"""
    if u in _cache: return _cache[u]
    v=list(base); v[u]=(v[u]+1)%P; forwardp(v)
    d={}
    for a in CHK:
        x=evalp(L.polys[a],v)
        if x!=bd[a]: d[a]=(x-bd[a])%P
    _cache[u]=d
    return d

def gf_solve(M, rhs, n):
    m=len(M)
    A=[list(M[i])+[rhs[i]] for i in range(m)]
    piv={}
    r=0
    for c in range(n):
        pr=next((i for i in range(r,m) if A[i][c]),None)
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]
        inv=pow(A[r][c],-1,P)
        A[r]=[x*inv%P for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                f=A[i][c]
                Ar=A[r]; Ai=A[i]
                A[i]=[(Ai[k]-f*Ar[k])%P for k in range(n+1)]
        piv[c]=r; r+=1
        if r==m: break
    for i in range(r,m):
        if A[i][n] and not any(A[i][:n]): return None
    x=[0]*n
    for c,i in piv.items(): x[c]=A[i][n]
    return x

def solve_closure(base, maxit=12, verbose=True):
    bd=allchk(base)
    F=failing(bd)
    if verbose: print(f"  start failing checks mod p: {len(F)} {F}")
    knobs=[]
    rows=set(F)
    for it in range(maxit):
        # knobs = free inputs touching any row
        cand=set()
        for a in rows:
            cand |= {u for u in L.avars[a] if u in FS}
        # also inputs one gate up: free inputs whose response touches a row
        newk=[u for u in sorted(cand) if u not in knobs]
        for u in newk:
            d=response(base,bd,u)
            if d: knobs.append(u); rows |= set(d)
        cols=[response(base,bd,u) for u in knobs]
        rowl=sorted(rows)
        idx={a:i for i,a in enumerate(rowl)}
        M=[[0]*len(knobs) for _ in rowl]
        for j,d in enumerate(cols):
            for a,x in d.items(): M[idx[a]][j]=x
        rhs=[(-bd[a])%P for a in rowl]
        t0=time.time()
        x=gf_solve(M,rhs,len(knobs))
        if verbose:
            print(f"  it{it}: system {len(rowl)}x{len(knobs)} -> "
                  f"{'SOLVED' if x is not None else 'inconsistent'} ({time.time()-t0:.0f}s)", flush=True)
        if x is None: return None,knobs,rowl
        v=list(base)
        for j,u in enumerate(knobs): v[u]=(v[u]+x[j])%P
        forwardp(v)
        d2=allchk(v); F2=failing(d2)
        if verbose: print(f"       applied -> failing checks mod p = {len(F2)} {F2[:14]}", flush=True)
        if not F2: return v,knobs,rowl
        # Newton: re-base and repeat
        base=v; bd=d2; _cache.clear(); rows=set(F2)|set(rows); knobs=[]
    return None,knobs,sorted(rows)

if __name__=='__main__':
    base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
    bd=allchk(base)
    # free wins first: x9118 kills a35759, x8731 kills a35760, neither disturbs anything
    for a,u in [(35759,9118),(35760,8731)]:
        x=solvep(a,u,base)
        if x is not None: base[u]=x
    forwardp(base)
    print("after the two free knobs:", failing(allchk(base)))
    t0=time.time()
    v,knobs,rows=solve_closure(base)
    print(f"total {time.time()-t0:.0f}s")
    if v is not None:
        json.dump([int(x) for x in v], open(os.path.join(HERE,'data','gmp5_state.json'),'w'))
        print("*** EVERY CHECK ZERO MOD P -- saved data/gmp5_state.json")
