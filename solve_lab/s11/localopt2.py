"""Exact local optimum, corrected.

A = atoms we allow to move (those appearing in the equations that touch the broken atoms).
KNOBS = variables occurring in no atom outside A -- moving them cannot touch anything else.
OBJ = every equation containing any atom of A (not just the ones touching the broken atoms;
that was the bug: an atom of A can appear in equations far outside).
Drop candidates are the currently-FAILING equations of OBJ, everything else must be preserved.
"""
import sys, os, json, time, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
from zsolve import solve_int
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)

def build(src, verbose=True):
    v=load_raw(src)
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    BR=[a for a in range(L.NA) if AV[a]!=0]
    E0=sorted(set().union(*[set(L.atom2eq.get(a,{})) for a in BR]))
    A=sorted(set().union(*[set(L.eq_atoms[e][2]) for e in E0]))
    AS=set(A)
    OBJ=sorted(set().union(*[set(L.atom2eq.get(a,{})) for a in A]))
    VARS=sorted(set().union(*[set(L.avars[a]) for a in A]))
    KN=[u for u in VARS if all(b in AS for b in L.var_atoms[u])]
    def inner(e,AVx): return sum(c*AVx[a] for a,c in L.eq_atoms[e][2].items())
    b=[inner(e,AV) for e in OBJ]
    FAILI=[i for i,x in enumerate(b) if x]
    if verbose:
        print(f"{os.path.basename(src)}: broken atoms={BR}")
        print(f"  movable atoms={len(A)}  objective equations={len(OBJ)}  failing={len(FAILI)}  knobs={len(KN)}")
    cols=[];used=[];nonaff=[]
    for u in KN:
        r=[]
        for step in (1,2):
            v[u]+=step
            A2=list(AV)
            for a in AS: A2[a]=L.evalpoly(L.polys[a],v)
            r.append([inner(e,A2)-b[i] for i,e in enumerate(OBJ)])
            v[u]-=step
        if all(r[1][i]==2*r[0][i] for i in range(len(OBJ))):
            if any(r[0]): cols.append(r[0]); used.append(u)
        else: nonaff.append(u)
    if verbose: print(f"  affine knobs with effect: {len(used)}  non-affine: {len(nonaff)}")
    M=[[cols[j][i] for j in range(len(cols))] for i in range(len(OBJ))]
    rhs=[-x for x in b]
    return v,OBJ,FAILI,used,M,rhs

if __name__=='__main__':
    src=sys.argv[1]
    v,OBJ,FAILI,used,M,rhs=build(src)
    t0=time.time()
    hit=None
    for nd in range(0,len(FAILI)+1):
        for D in itertools.combinations(FAILI,nd):
            keep=[i for i in range(len(OBJ)) if i not in set(D)]
            x=solve_int([M[i] for i in keep],[rhs[i] for i in keep])
            if x is not None: hit=(nd,D,x); break
        print(f"    drop {nd}: {'FEASIBLE' if hit else 'infeasible'} ({time.time()-t0:.0f}s)", flush=True)
        if hit: break
    if hit:
        nd,D,x=hit
        print(f"  MINIMUM local failures = {nd}  drop {[OBJ[i] for i in D]}")
        for j,u in enumerate(used): v[u]+=x[j]
        AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
        F=L.failing_eqs(AV)
        print(f"APPLIED -> failing={len(F)} score={L.NEQ-len(F)}  broken atoms={[a for a in range(L.NA) if AV[a]!=0]}")
        if len(F)<7:
            out=os.path.join(HERE,'data','localopt2_out.json')
            json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(out,'w')); print("saved",out)
