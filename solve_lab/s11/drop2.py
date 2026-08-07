"""Which rows of the ripple-response system carry the p-obstruction, and what do they cost?

The system is rationally consistent but not integrally (nor mod p) solvable, so a small set of
check atoms must stay broken.  The SCORE cost of leaving atom a broken is the number of
equations containing a -- and many of the big aggregate checks sit in exactly one equation.
So: search for the cheapest drop set, ordered by cost rather than by cardinality.
"""
import sys, os, json, time, itertools, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
from zsolve import solve_int
import resp as R
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
src=sys.argv[1]; hops=int(sys.argv[2]) if len(sys.argv)>2 else 2
v=load_raw(src)
BR=[a for a in range(L.NA) if R.av(v,a)!=0]
C=R.candidates(v,BR,hops)
cols={}
for u in C:
    d1,_,_=R.response(v,u,1)
    if not d1: continue
    d2,_,_=R.response(v,u,2)
    if all(d2.get(a,0)==2*d1.get(a,0) for a in set(d1)|set(d2)): cols[u]=d1
ROWS=sorted(set(BR)&R.ISCHK | set().union(*[set(d) for d in cols.values()]))
used=sorted(cols)
M=[[cols[u].get(a,0) for u in used] for a in ROWS]
rhs=[-R.av(v,a) for a in ROWS]
NEQ=[len(L.atom2eq.get(a,{})) for a in ROWS]
print(f"system {len(ROWS)}x{len(used)}   rows by #equations: {dict(sorted(collections.Counter(NEQ).items())[:8])}")
cheap=sorted(range(len(ROWS)), key=lambda i:(NEQ[i],i))
print("cheapest rows:", [(ROWS[i],NEQ[i]) for i in cheap[:20]])
t0=time.time()
def feas(drop):
    keep=[i for i in range(len(ROWS)) if i not in drop]
    return solve_int([M[i] for i in keep],[rhs[i] for i in keep])
best=None
for k in (1,2,3):
    pool=cheap[:40] if k>1 else cheap
    found=False
    for D in itertools.combinations(pool,k):
        cost=len(set().union(*[set(L.atom2eq.get(ROWS[i],{})) for i in D]))
        if best is not None and cost>=best[0]: continue
        x=feas(set(D))
        if x is not None:
            best=(cost,D,x); found=True
            print(f"  drop {[ROWS[i] for i in D]} cost={cost} eqs  ({time.time()-t0:.0f}s)", flush=True)
    print(f"  finished k={k} best={best[0] if best else None} ({time.time()-t0:.0f}s)", flush=True)
    if best and best[0]<=2: break
if best:
    cost,D,x=best
    print(f"BEST drop set {[ROWS[i] for i in D]} -> {cost} failing equations from these rows")
    seeds={u:v[u]+x[j] for j,u in enumerate(used) if x[j]}
    L.ripple(v,seeds)
    AV=[R.av(v,a) for a in range(L.NA)]
    F=L.failing_eqs(AV)
    print(f"APPLIED -> broken atoms={len([a for a in range(L.NA) if AV[a]!=0])} failing={len(F)} score={L.NEQ-len(F)}")
    json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','drop2_out.json'),'w'))
