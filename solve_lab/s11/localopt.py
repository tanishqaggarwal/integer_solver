"""Generic exact local optimum, the local1 computation for any state.

Take the broken atoms, the equations touching them, and the variables that occur in NO atom
outside those equations (fully local knobs, so moving them cannot disturb anything else).
Measure each knob's exact effect on every equation's inner value, verify it is affine, and then
enumerate drop-sets by size: the smallest feasible one is the exact local minimum of failing
equations.
"""
import sys, os, json, time, itertools, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
from zsolve import solve_int
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)

def run(src, extra=(), verbose=True):
    v=load_raw(src)
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    BR=[a for a in range(L.NA) if AV[a]!=0]
    EQS=sorted(set().union(*[set(L.atom2eq.get(a,{})) for a in BR]))
    ATOMS=sorted(set().union(*[set(L.eq_atoms[e][2]) for e in EQS]))
    AS=set(ATOMS)
    # every equation containing ANY atom we allow to move is part of the objective
    EQS=sorted(set().union(*[set(L.atom2eq.get(a,{})) for a in ATOMS]))
    VARS=sorted(set().union(*[set(L.avars[a]) for a in ATOMS]))
    LOCAL=[u for u in VARS if all(a in AS for a in L.var_atoms[u])]
    LOCAL=sorted(set(LOCAL)|set(extra))
    def inner(e, AV): return sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())
    b=[inner(e,AV) for e in EQS]
    if verbose:
        print(f"{os.path.basename(src)}: broken={len(BR)} {BR}  equations={len(EQS)}  "
              f"atoms={len(ATOMS)}  local knobs={len(LOCAL)}")
        print(f"  currently failing: {sum(1 for x in b if x)} of {len(EQS)}")
    cols=[]; used=[]; nonaff=[]
    for u in LOCAL:
        r=[]
        ok=True
        for step in (1,2):
            v[u]+=step
            A2=[L.evalpoly(L.polys[a],v) if a in AS else 0 for a in range(L.NA)]
            for a in range(L.NA):
                if a not in AS: A2[a]=AV[a]
            r.append([inner(e,A2)-b[i] for i,e in enumerate(EQS)])
            v[u]-=step
        if all(r[1][i]==2*r[0][i] for i in range(len(EQS))):
            if any(r[0]): cols.append(r[0]); used.append(u)
        else:
            nonaff.append(u)
    if verbose: print(f"  affine knobs with an effect: {len(used)}   non-affine: {len(nonaff)} {nonaff[:10]}")
    M=[[cols[j][i] for j in range(len(cols))] for i in range(len(EQS))]
    rhs=[-x for x in b]
    t0=time.time()
    for nd in range(0, len(EQS)+1):
        for D in itertools.combinations(range(len(EQS)), nd):
            keep=[i for i in range(len(EQS)) if i not in D]
            x=solve_int([M[i] for i in keep],[rhs[i] for i in keep])
            if x is not None:
                print(f"  MINIMUM local failures = {nd}   drop {[EQS[i] for i in D]}  ({time.time()-t0:.0f}s)")
                return v, EQS, used, x, D
        if verbose: print(f"    drop {nd}: infeasible ({time.time()-t0:.0f}s)", flush=True)
    return v, EQS, used, None, None

if __name__=='__main__':
    src=sys.argv[1]
    v,EQS,used,x,D=run(src)
    if x is not None:
        for j,u in enumerate(used): v[u]+=x[j]
        AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
        F=L.failing_eqs(AV)
        print(f"APPLIED -> failing={len(F)} score={L.NEQ-len(F)}")
        if len(F)<7:
            out=os.path.join(HERE,'data','localopt_out.json')
            json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(out,'w')); print("saved",out)
