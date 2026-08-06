"""Minimum-cost achievable residual, mod p.

The non-bit system is inconsistent, so some checks must stay nonzero.  But WHICH ones is our
choice, and the cost of a check is the number of EQUATIONS containing it -- and several of the
checks in the obstruction certificates sit in exactly one equation.  So: drop candidate rows
(allow them to be nonzero), test consistency of the rest, and minimise total equation cost.
"""
import sys, os, json, time, pickle, itertools, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
from gmp7 import failing
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
D=pickle.load(open(os.path.join(HERE,'data','resp_modp.pkl'),'rb'))
base=D['base']; bd=D['bd']; cols=D['cols']
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
knobs=sorted(u for u in cols if u not in BITS)
F=failing(bd)
rows=set(F)
for u in knobs: rows |= set(cols[u])
rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
COST={a:len(L.atom2eq.get(a,{})) for a in rowl}
n=len(knobs)
COL=[{idx[a]:x%P for a,x in cols[u].items()} for u in knobs]
RHS={i:(-bd[rowl[i]])%P for i in range(len(rowl)) if bd.get(rowl[i])}
def consistent(drop):
    """is the system solvable when rows in `drop` are unconstrained?"""
    keep=[i for i in range(len(rowl)) if i not in drop]
    kidx={i:t for t,i in enumerate(keep)}
    m=len(keep)
    R=[dict() for _ in range(m)]
    for j,c in enumerate(COL):
        for i,x in c.items():
            if i in kidx: R[kidx[i]][j]=x
    b=[RHS.get(i,0) for i in keep]
    used=[False]*m
    for c in range(n):
        cand=[i for i in range(m) if not used[i] and c in R[i]]
        if not cand: continue
        i=min(cand,key=lambda i:len(R[i])); used[i]=True
        inv=pow(R[i][c],-1,P)
        R[i]={k:x*inv%P for k,x in R[i].items()}; b[i]=b[i]*inv%P
        Ri=R[i]; bi=b[i]
        for k in cand:
            if k==i: continue
            f=R[k].get(c)
            if not f: continue
            Rk=R[k]
            for kk,x in Ri.items():
                nv=(Rk.get(kk,0)-f*x)%P
                if nv: Rk[kk]=nv
                elif kk in Rk: del Rk[kk]
            b[k]=(b[k]-f*bi)%P
    return not any((not R[i]) and b[i] for i in range(m))
print(f"system {len(rowl)} rows x {n} knobs; failing {F}")
print("no drops consistent:", consistent(set()))
CERT=set([2423,3576,3578,7930,19297,19299,21617,25676,29539,30984,31670,33792,40562,40623,
          26731,31672,33929,36185,40812,42245])
pool=sorted(CERT, key=lambda a:(COST[a],a))
print("certificate rows by equation cost:", [(a,COST[a]) for a in pool])
t0=time.time(); best=None
for k in (1,2,3,4):
    for S in itertools.combinations(pool,k):
        cost=len(set().union(*[set(L.atom2eq.get(a,{})) for a in S]))
        if best is not None and cost>=best[0]: continue
        if consistent({idx[a] for a in S if a in idx}):
            best=(cost,S)
            print(f"  FEASIBLE dropping {S}  -> {cost} failing equations  ({time.time()-t0:.0f}s)", flush=True)
    print(f"  finished k={k} best={best[0] if best else None} ({time.time()-t0:.0f}s)", flush=True)
    if best and best[0]<=2: break
print("BEST:",best)
if best: json.dump({'cost':best[0],'drop':list(best[1])}, open(os.path.join(HERE,'data','gmp25.json'),'w'))
