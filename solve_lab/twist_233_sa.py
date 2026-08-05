#!/usr/bin/env python3
"""SA repair of the 233-SIDE twist activation (freeze x_18274:=x_9770, x_17728:=
x_3183). This core is SQUARE-FREE (verifier squares hold automatically since the
twist holds); the ~12 broken atoms are the x_8821 denominator web (a4954:
x_8821*x_18274=x_6773, a13204, a41285, ...). Different, possibly more tractable
frustrated core than the slack-active one. Min-conflicts + SA over integer vars."""
import json, time, random, math, sys
from collections import defaultdict
from confluent_eval5 import build5, make_forward
from twist_233side import make_frz_solver
from slack_active import viol_atoms
from propagate import atom_vars, NVARS

def resid(poly, val):
    s = 0
    for m, c in poly.items():
        t = c
        for x in m: t *= val[x]
        s += t
    return s

def solve_for(poly, val, v):
    c0=c1=c2=0
    for m,c in poly.items():
        k=m.count(v); t=c
        for x in m:
            if x!=v: t*=val[x]
        if k==0: c0+=t
        elif k==1: c1+=t
        else: c2+=t
    if c2==0:
        if c1==0 or (-c0)%c1: return None
        return (-c0)//c1
    if c1==0:
        if (-c0)%c2: return None
        r=(-c0)//c2
        if r<0: return None
        s=math.isqrt(r); return s if s*s==r else None
    disc=c1*c1-4*c2*c0
    if disc<0: return None
    s=math.isqrt(disc)
    if s*s!=disc: return None
    for num in (-c1+s,-c1-s):
        if num%(2*c2)==0: return num//(2*c2)
    return None

def main():
    t0=time.time()
    A,kind,info,seq0,bestval,ncyc=build5()
    order=json.load(open('eval_order.json'))['order']
    defset=set(v for v in kind if kind[v]!='const')
    seq=[v for v in order if v in defset and v not in (9770,3183)]
    seq+=[v for v in (9770,3183) if v in defset]
    seq+=[v for v in defset if v not in set(order) and v not in (9770,3183)]
    solve=make_forward(kind,info,seq,bestval)
    run=make_frz_solver(kind,info,seq,bestval,{18274,17728})
    ACT=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 and sys.argv[1] else []
    SEED=int(sys.argv[2]) if len(sys.argv)>2 else 7
    OUT=sys.argv[3] if len(sys.argv)>3 else 'cand_233_SOLVED.json'

    var_atoms=defaultdict(list)
    for a,poly in enumerate(A):
        for v in atom_vars(poly): var_atoms[v].append(a)

    v1=solve(list(bestval),ACT)
    frozen={18274:v1[9770], 17728:v1[3183]}
    val=run(list(v1),frozen)
    bad=set(a for a in range(len(A)) if resid(A[a],val))
    print(f"seed act={ACT}: {len(bad)} violated: {sorted(bad)} ({time.time()-t0:.0f}s)", flush=True)
    PROT={9770,3183,18274,17728}
    rng=random.Random(SEED)
    best_bad=len(bad); T=3.0; it=0; TIME=5400
    while time.time()-t0<TIME and best_bad>0:
        it+=1
        if not bad: break
        a=rng.choice(tuple(bad)); poly=A[a]
        cand=[v for v in atom_vars(poly) if v not in PROT]
        if not cand: bad.discard(a); continue
        v=rng.choice(cand); nv=solve_for(poly,val,v)
        if nv is None or nv==val[v]: continue
        old=val[v]
        before=sum(1 for aa in var_atoms[v] if resid(A[aa],val))
        val[v]=nv
        after=sum(1 for aa in var_atoms[v] if resid(A[aa],val))
        d=after-before
        if d<=0 or rng.random()<math.exp(-d/max(T,0.05)):
            for aa in var_atoms[v]:
                if resid(A[aa],val): bad.add(aa)
                else: bad.discard(aa)
        else:
            val[v]=old
        if len(bad)<best_bad:
            best_bad=len(bad)
            print(f"  it {it}: NEW BEST {best_bad}: {sorted(bad)} (T={T:.2f}, {time.time()-t0:.0f}s)", flush=True)
            if best_bad==0:
                allbad=viol_atoms(A,val)
                print(f"  full verify: {len(allbad)}", flush=True)
                if not allbad:
                    json.dump({f"x_{i}":val[i] for i in range(NVARS)},open(OUT,'w'))
                    print(f"  *** SOLVED *** -> {OUT}", flush=True); return
                bad=set(allbad); best_bad=len(bad)
        if it%4000==0:
            T*=0.9
            if T<0.1: T=3.0
        if it%50000==0:
            print(f"  it {it}: bad={len(bad)} best={best_bad} T={T:.2f} ({time.time()-t0:.0f}s)", flush=True)
    print(f"done: best {best_bad}: {sorted(bad)} ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
