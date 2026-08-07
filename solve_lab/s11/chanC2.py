"""Cheapest set of checks to leave nonzero in channel C.

Its certificates are small: {a688,a40608} has just two rows, and several members sit in ONE
equation.  Cost of a drop set = |union of the equations containing those atoms|, which is exactly
the score penalty, so search it by cost against the 7 the checkpoint achieves.
"""
import sys, os, json, time, pickle, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
P=L.P; sys.set_int_max_str_digits(400000)
D=pickle.load(open(os.path.join(HERE,'data','resp_C.pkl'),'rb'))
base=D['base']; bd=D['bd']; cols=D['cols']
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
knobs=sorted(u for u in cols if u not in BITS)
F=sorted(a for a,x in bd.items() if x)
rows=set(F)
for u in knobs: rows |= set(cols[u])
COLS=[cols[u] for u in knobs]
EQ=lambda a: set(L.atom2eq.get(a,{}))
def consistent(drop):
    rl=sorted(rows-set(drop)); idx={a:i for i,a in enumerate(rl)}
    m=len(rl); n=len(COLS)
    R=[dict() for _ in rl]
    for j,d in enumerate(COLS):
        for a,x in d.items():
            if a in idx: R[idx[a]][j]=x%P
    b=[(-bd[a])%P for a in rl]
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
CERTROWS=set([1618,2423,31670,2421,10506,688,26731,31672,33929,33937,33938,40068,40608,
              37640,40969,37638,42117,25676,33792,40562,40623,42245])
print("certificate rows and their equation counts:",
      sorted(((a,len(EQ(a))) for a in CERTROWS), key=lambda z:z[1]))
pool=sorted(CERTROWS, key=lambda a: len(EQ(a)))
t0=time.time(); best=(7,None)
for k in range(1,6):
    for S in itertools.combinations(pool,k):
        cost=len(set().union(*[EQ(a) for a in S]))
        if cost>=best[0]: continue
        if consistent(set(S)):
            best=(cost,S); print(f"  *** cost {cost}: leave {S} nonzero ({time.time()-t0:.0f}s)", flush=True)
    print(f"  k={k} done best={best[0]} ({time.time()-t0:.0f}s)", flush=True)
print("BEST:",best)
