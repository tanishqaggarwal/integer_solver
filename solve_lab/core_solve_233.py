#!/usr/bin/env python3
"""Exact-ish solve of the 233-side residual core from a checkpoint state. Strategy:
prioritize FREE-variable moves (df=None vars appear in few atoms -> satisfy an atom
with minimal ripple), then least-ripple var moves, with limited backtracking. Loads
cand_233r_*_ckpt.json (val + bad), tries to drive ALL atoms to 0."""
import json, time, sys, math
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS

CK = sys.argv[1] if len(sys.argv) > 1 else 'cand_233r_55_ckpt.json'

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
    A=load_atoms()
    prov=json.load(open('eval_order.json'))['prov']
    freevars=set(v for v in range(NVARS) if not (v<len(prov) and prov[v] and prov[v][0]>=0))
    d=json.load(open(CK)); val=[0]*NVARS
    for k,x in d['val'].items(): val[int(k)]=x
    va=defaultdict(list)
    for a,poly in enumerate(A):
        for v in atom_vars(poly): va[v].append(a)
    def bad_of(v): return [aa for aa in va[v] if resid(A[aa],val)]
    bad=set(a for a in range(len(A)) if resid(A[a],val))
    print(f"loaded {CK}: {len(bad)} violated: {sorted(bad)} ({time.time()-t0:.0f}s)", flush=True)

    # product-def map: free var -> (out_var, product_atom) where out = coef*prod(factors)
    prod_def={}
    for a,poly in enumerate(A):
        if len(poly)==2:
            ms=sorted(poly.keys(),key=len)
            if len(ms[0])==1 and len(ms[1])==2:  # out - factorA*factorB
                out=ms[0][0]; f1,f2=ms[1]
                for fv,ov in [(f1,f2),(f2,f1)] if f1!=f2 else [(f1,f1)]:
                    if fv in freevars: prod_def.setdefault(fv,(out,a,ov))
    def joint_solve(a, fv):
        """solve broken atom a for free factor fv, with its product-def substituted."""
        out,pa,ov=prod_def[fv]; ppoly=A[pa]
        cout=ppoly[(out,)]; cprod=ppoly[tuple(sorted((fv,ov)))]  # cout*out + cprod*fv*ov=0 => out=-cprod*fv*ov/cout
        # build atom a as poly in fv after substituting out = k*fv where k=-cprod*ov/cout (needs cout|...)
        k_num=-cprod*val[ov]; k_den=cout
        if k_den==0 or k_num%k_den: return None
        kf=k_num//k_den   # out = kf*fv
        c0=c1=c2=0
        for m,c in A[a].items():
            t=c; deg=0
            for x in m:
                if x==fv: deg+=1
                elif x==out: deg+=1; t*=kf
                else: t*=val[x]
            if deg==0: c0+=t
            elif deg==1: c1+=t
            else: c2+=t
        if c2==0:
            if c1==0 or (-c0)%c1: return None
            return (-c0)//c1, out, kf
        disc=c1*c1-4*c2*c0
        if disc<0: return None
        s=math.isqrt(disc)
        if s*s!=disc: return None
        for num in (-c1+s,-c1-s):
            if num%(2*c2)==0: return num//(2*c2), out, kf
        return None

    import random
    rng=random.Random(0)
    best=len(bad); best_val=list(val); stall=0
    for it in range(2000000):
        if not bad: break
        a=rng.choice(tuple(bad)); poly=A[a]
        # TRY joint solve: free factor var in this atom
        jvars=[v for v in atom_vars(poly) if v in prod_def]
        if jvars and rng.random()<0.5:
            fv=rng.choice(jvars); r=joint_solve(a,fv)
            if r is not None:
                nvfv,out,kf=r; oldfv=val[fv]; oldout=val[out]
                val[fv]=nvfv; val[out]=kf*nvfv
                touched=set(va[fv])|set(va[out])
                b0=None  # accept if it doesn't increase local viol
                for aa in touched:
                    if resid(A[aa],val): bad.add(aa)
                    else: bad.discard(aa)
                continue
        cvars=[v for v in atom_vars(poly) if v not in {9770,3183,18274,17728}]
        # prefer free vars, then fewest-atom vars
        cvars.sort(key=lambda v: (v not in freevars, len(va[v])))
        moved=False
        for v in cvars[:8]:
            nv=solve_for(poly,val,v)
            if nv is None or nv==val[v]: continue
            old=val[v]
            b0=sum(1 for aa in va[v] if resid(A[aa],val))
            val[v]=nv
            b1=sum(1 for aa in va[v] if resid(A[aa],val))
            if b1<=b0 or rng.random()<0.25*math.exp(-(b1-b0)):
                for aa in va[v]:
                    if resid(A[aa],val): bad.add(aa)
                    else: bad.discard(aa)
                moved=True; break
            else:
                val[v]=old
        if not moved:
            v=rng.choice(cvars); nv=solve_for(poly,val,v)
            if nv is not None:
                val[v]=nv
                for aa in va[v]:
                    if resid(A[aa],val): bad.add(aa)
                    else: bad.discard(aa)
        if len(bad)<best:
            best=len(bad); best_val=list(val); stall=0
            print(f"  it {it}: NEW BEST {best}: {sorted(bad)} ({time.time()-t0:.0f}s)", flush=True)
            if best==0:
                allbad=[a for a in range(len(A)) if resid(A[a],val)]
                print(f"  verify: {len(allbad)}", flush=True)
                if not allbad:
                    json.dump({f"x_{i}":val[i] for i in range(NVARS)},open('cand_CORESOLVED.json','w'))
                    print("  *** SOLVED ***", flush=True); return
        else:
            stall+=1
            if stall>50000:
                val=list(best_val); bad=set(a for a in range(len(A)) if resid(A[a],val)); stall=0
        if it%200000==0 and it:
            print(f"  it {it}: bad={len(bad)} best={best} ({time.time()-t0:.0f}s)", flush=True)
    print(f"done: best {best} ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
