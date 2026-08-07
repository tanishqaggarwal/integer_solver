"""EXACT ripple-response system.

Key simplification learned today: a gate ripple keeps every GATE atom at zero automatically,
so only CHECK atoms need to appear as rows.  For each candidate FREE variable u we measure the
exact effect of u -> u+1 and u -> u+2 after a ripple; if the second is twice the first the
response is exactly affine and the column is trusted (we record, not silently drop, the rest).

Then solve   R . delta  =  -(current check-atom values)   over Z.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
from zsolve import solve_int
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
CHECK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
ISCHK=set(CHECK)

def av(v,a): return L.evalpoly(L.polys[a],v)

def response(v, u, step):
    v2=list(v); L.ripple(v2,{u:v[u]+step})
    d={}
    diff=[t for t in range(L.NVARS) if v2[t]!=v[t]]
    cand=set()
    for t in diff: cand|=set(L.var_atoms[t])
    for a in cand:
        if a in ISCHK:
            n=av(v2,a)-av(v,a)
            if n: d[a]=n
    broke_gate=[a for a in cand if a not in ISCHK and av(v2,a)!=0]
    return d, broke_gate, len(diff)

def candidates(v, seedatoms, hops=2):
    C=set(); frontier=set(seedatoms)
    seen=set()
    for h in range(hops):
        nxt=set()
        for a in frontier:
            if a in seen: continue
            seen.add(a)
            for t in L.avars[a]:
                if t in FREE: C.add(t)
                for b in L.var_atoms[t]: nxt.add(b)
        frontier=nxt
    return sorted(C)

if __name__=='__main__':
    src=sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
    hops=int(sys.argv[2]) if len(sys.argv)>2 else 1
    v=load_raw(src)
    BR=[a for a in range(L.NA) if av(v,a)!=0]
    print(f"{os.path.basename(src)}: broken atoms {BR}")
    C=candidates(v, BR, hops)
    print(f"candidates ({hops} hops): {len(C)} free variables")
    t0=time.time()
    cols={}; nonlin=[]; gatebreak=[]
    for u in C:
        d1,gb1,n1=response(v,u,1)
        if not d1 and not gb1: continue
        d2,gb2,n2=response(v,u,2)
        ok=all(d2.get(a,0)==2*d1.get(a,0) for a in set(d1)|set(d2))
        if gb1 or gb2: gatebreak.append(u)
        if not ok: nonlin.append(u); continue
        cols[u]=d1
    print(f"  affine columns: {len(cols)}  non-affine: {len(nonlin)} {nonlin[:15]}  "
          f"ripple-unrepairable gates from: {len(gatebreak)} {gatebreak[:10]}  ({time.time()-t0:.0f}s)")
    ROWS=sorted(set(BR)&ISCHK | set().union(*[set(d) for d in cols.values()]) if cols else set(BR)&ISCHK)
    used=sorted(cols)
    M=[[cols[u].get(a,0) for u in used] for a in ROWS]
    rhs=[-av(v,a) for a in ROWS]
    print(f"  system {len(ROWS)} check rows x {len(used)} cols")
    t1=time.time()
    x=solve_int(M,rhs)
    print(f"  integer-solvable: {'YES' if x is not None else 'no'} ({time.time()-t1:.0f}s)")
    if x is not None:
        seeds={u:v[u]+x[j] for j,u in enumerate(used) if x[j]}
        print("  seeds:",{k:('%+d'%s if abs(s)<10**8 else 'big') for k,s in [(u,x[j]) for j,u in enumerate(used) if x[j]]})
        L.ripple(v,seeds)
        AV=[av(v,a) for a in range(L.NA)]
        F=L.failing_eqs(AV)
        print(f"  APPLIED -> broken atoms={len([a for a in range(L.NA) if AV[a]!=0])} failing={len(F)} score={L.NEQ-len(F)}")
        json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','resp_out.json'),'w'))
