"""General repair loop.

After the neighbourhood fix, every surviving broken atom turned out to be a CHECK of the form
"a free variable must equal a computed value".  So: repeatedly pick a broken atom, find a
variable that occurs linearly in it and whose adjustment zeroes it, prefer the one that
disturbs the fewest other atoms, set it, and ripple the gates.  Repeat to a fixpoint.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)

def lin_solve(a, t, v):
    """value of t zeroing atom a (t must occur linearly), else None"""
    c=0
    for m,cc in L.polys[a].items():
        k=m.count(t)
        if k==0: continue
        if k>1: return None
        term=cc
        for u in m:
            if u!=t: term*=v[u]
        c+=term
    if c==0: return None
    old=v[t]; v[t]=0
    rest=L.evalpoly(L.polys[a], v)
    v[t]=old
    if rest % c: return None
    return -rest//c

def nonzero(v):
    return [a for a in range(L.NA) if L.evalpoly(L.polys[a],v)!=0]

def repair(v, rounds=60, verbose=True):
    for it in range(rounds):
        NZ=nonzero(v)
        if not NZ:
            return v, []
        best=None
        for a in NZ:
            for t in sorted(L.avars[a], key=lambda t: len(L.var_atoms[t])):
                x=lin_solve(a,t,v)
                if x is None or x==v[t]: continue
                cost=len(L.var_atoms[t])-1
                if best is None or cost<best[0]:
                    best=(cost,a,t,x)
                if cost==0: break
            if best and best[0]==0: break
        if best is None:
            if verbose: print(f"   it{it}: stuck, {len(NZ)} broken atoms {NZ[:12]}")
            return v, NZ
        cost,a,t,x=best
        L.ripple(v, {t:x})
        if verbose:
            NZ2=nonzero(v)
            print(f"   it{it}: fixed a{a} via x{t} (cost {cost}) -> broken atoms {len(NZ2)} {NZ2[:12]}", flush=True)
    return v, nonzero(v)

if __name__=='__main__':
    src=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'data','fix2_round.json')
    v=load_raw(src)
    NZ=nonzero(v)
    AV=[atomval(a,v) for a in range(L.NA)]
    print(f"start: broken atoms={len(NZ)} {NZ}  failing eqs={len(L.failing_eqs(AV))}")
    t0=time.time()
    v,NZ=repair(v)
    AV=[atomval(a,v) for a in range(L.NA)]
    F=L.failing_eqs(AV)
    print(f"RESULT: broken atoms={len(NZ)} failing eqs={len(F)} score={L.NEQ-len(F)}  ({time.time()-t0:.0f}s)")
    out=os.path.join(HERE,'data','repair_out.json')
    json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(out,'w'))
    print("saved", out)
