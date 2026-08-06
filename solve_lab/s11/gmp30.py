"""Cheapest combination: break some cheap gates AND leave some cheap checks nonzero.

Score penalty = |equations containing the broken gate atoms UNION the equations containing the
checks left nonzero|.  The 39,026 checkpoint is one point in this space (7 equations); the
question is whether a better one exists.  Each consistency test reuses the cached columns, so
thousands of combinations are affordable.
"""
import sys, os, json, time, pickle, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
Z=pickle.load(open(os.path.join(HERE,'data','gmp29_cols.pkl'),'rb'))
C=Z['C']; G=Z['G']
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
bd={a:evalp(L.polys[a],base) for a in CHK}
F=[a for a in CHK if bd[a]]
knobcols=[C[u] for u in sorted(C)]
EQ=lambda a: set(L.atom2eq.get(a,{}))
def consistent(gates, drop):
    cols=knobcols+[G[g] for g in gates]
    rows=set(F)
    for d in cols: rows |= set(d)
    rows-=set(drop)
    rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
    m=len(rowl); n=len(cols)
    R=[dict() for _ in rowl]
    for j,d in enumerate(cols):
        for a,x in d.items():
            if a in idx: R[idx[a]][j]=x
    b=[(-bd[a])%P for a in rowl]
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
GC=sorted(G, key=lambda g: len(EQ(g)))
allrows=set(F)
for d in knobcols: allrows |= set(d)
for d in G.values(): allrows |= set(d)
cheapchk=sorted(allrows, key=lambda a: len(EQ(a)))
print(f"gate candidates {len(GC)} (cheapest {[(g,len(EQ(g))) for g in GC[:6]]})")
print(f"check candidates {len(cheapchk)} (cheapest {[(a,len(EQ(a))) for a in cheapchk[:8]]})")
t0=time.time(); best=(7,None,None)     # must beat the checkpoint's 7
tested=0
for ng in range(0,3):
    for gs in itertools.combinations(GC[:14], ng):
        gcost=set().union(*[EQ(g) for g in gs]) if gs else set()
        if len(gcost)>=best[0]: continue
        for nd in range(0,4):
            for ds in itertools.combinations(cheapchk[:16], nd):
                cost=len(gcost | (set().union(*[EQ(a) for a in ds]) if ds else set()))
                if cost>=best[0]: continue
                tested+=1
                if consistent(gs, set(ds)):
                    best=(cost,gs,ds)
                    print(f"  *** cost {cost}: break {gs}, leave {ds} nonzero ({time.time()-t0:.0f}s)", flush=True)
    print(f"  ng={ng} done, best={best[0]}, tested={tested} ({time.time()-t0:.0f}s)", flush=True)
print("BEST:",best)
json.dump({'cost':best[0],'gates':list(best[1] or []),'drop':list(best[2] or [])},
          open(os.path.join(HERE,'data','gmp30.json'),'w'))
