"""Direct mod-p closure with upstream-cone knob discovery.

x9118 -> a35759 alone and x8731 -> a35760 alone, so two of the six fall for free.
x14853 -> {a29539, a40826} + {a19297, a19299, a30984, a36185, a40812}
x24548 -> {a7930,  a41512} + {a21617, a37662}
The collateral atoms carry no free variable themselves, but their inputs do, a gate or two
upstream.  So discover knobs by walking the definition DAG upward from each row atom, and solve
the whole thing at once over GF(p).
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
from gmp7 import allchk, failing, resp_at, sparse_solve
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)

def cone_free(atoms, depth):
    """free inputs reachable by walking upward from these atoms' variables"""
    seen=set(); frontier=set()
    for a in atoms: frontier |= set(L.avars[a])
    out=set()
    for _ in range(depth):
        nxt=set()
        for u in frontier:
            if u in seen: continue
            seen.add(u)
            if u in FREE: out.add(u); continue
            d=L.definer.get(u)
            if d is None: continue
            nxt |= set(L.avars[d])
        if not nxt: break
        frontier=nxt
    return out

def solve(base, depth=4, maxnewton=6, verbose=True):
    for nt in range(maxnewton):
        bd=allchk(base); F=failing(bd)
        if verbose: print(f"NEWTON {nt}: failing = {len(F)} {F[:16]}", flush=True)
        if not F: return base
        rows=set(F); knobs=set()
        for g in range(8):
            nk=cone_free(rows, depth) - knobs
            if not nk and knobs: break
            knobs |= nk
            R2=set(rows)
            for u in sorted(knobs):
                d=resp_at(base,bd,u)
                R2 |= set(d)
            if verbose: print(f"   grow {g}: knobs={len(knobs)} rows={len(R2)}", flush=True)
            if R2==rows: break
            rows=R2
            if len(rows)>1200 or len(knobs)>1200: break
        rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
        knobs=sorted(knobs)
        t0=time.time(); cols=[]
        for u in knobs:
            d=resp_at(base,bd,u,rowl)
            cols.append({idx[a]:x for a,x in d.items()})
        rhs=[(-bd[a])%P for a in rowl]
        if verbose: print(f"   responses {time.time()-t0:.0f}s; solving {len(rowl)}x{len(knobs)}", flush=True)
        t0=time.time()
        x=sparse_solve(cols,rowl,rhs)
        if verbose: print(f"   -> {'SOLVED' if x is not None else 'INCONSISTENT'} ({time.time()-t0:.0f}s)", flush=True)
        if x is None: return None
        v=list(base)
        for j,u in enumerate(knobs): v[u]=(v[u]+x[j])%P
        forwardp(v); base=v
    return None

if __name__=='__main__':
    base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
    depth=int(sys.argv[1]) if len(sys.argv)>1 else 4
    t0=time.time()
    v=solve(base, depth=depth)
    print(f"total {time.time()-t0:.0f}s")
    if v is not None:
        json.dump([int(x) for x in v], open(os.path.join(HERE,'data','gmp9_state.json'),'w'))
        print("*** EVERY CHECK ZERO MOD P -- saved data/gmp9_state.json")
