"""Targeted: the cheap checks among the obstruction are a40826 and a41512 -- ONE equation each.

If a7930 and a29539 can be zeroed once those two are released, the cost is 2 equations.
Test the promising drop sets directly instead of enumerating five million clusters.
"""
import sys, os, json, time, pickle, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
from gmp26 import forwardp_frozen
P=L.P; sys.set_int_max_str_digits(400000)
CHKS=set(a for a in range(L.NA) if L.atom_out.get(a) is None)
EQ=lambda a: set(L.atom2eq.get(a,{}))
Z=pickle.load(open(os.path.join(HERE,'data','gmp29_cols.pkl'),'rb')); C=Z['C']; G=Z['G']
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
bd={a:evalp(L.polys[a],base) for a in CHKS}
F=[a for a in CHKS if bd[a]]
print("failing checks:",[(a,len(EQ(a))) for a in F])
knobcols=[C[u] for u in sorted(C)]
def consistent(cols, drop):
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
def colfor(t, frozen):
    v=list(base); v[t]=(v[t]+1)%P; forwardp_frozen(v,frozen)
    d={a:(evalp(L.polys[a],v)-bd[a])%P for a in CHKS}
    return {a:x for a,x in d.items() if x}
ONEQ=sorted([a for a in CHKS if len(EQ(a))==1], key=lambda a:a)
print("1-equation checks:",len(ONEQ))
CANDDROP=[40826,41512,25676,42245,36185,40812,37662,40623,40562,33792]
GCAND=[36244,36245,36246,34869]
tests=[]
for nd in range(1,6):
    for ds in itertools.combinations(CANDDROP,nd):
        tests.append((tuple(),ds))
for ng in (1,2):
    for gs in itertools.combinations(GCAND,ng):
        for nd in range(0,5):
            for ds in itertools.combinations(CANDDROP,nd):
                tests.append((gs,ds))
def cost(gs,ds):
    s=set()
    for g in gs: s|=EQ(g)
    for a in ds: s|=EQ(a)
    return len(s)
tests=[t for t in tests if cost(*t)<7]
tests.sort(key=lambda t: cost(*t))
print(f"{len(tests)} configurations with cost < 7")
t0=time.time(); best=None
for gs,ds in tests:
    frozen={L.atom_out[g][1] for g in gs}
    cols=knobcols+[colfor(t,frozen) for t in frozen]
    if consistent(cols,set(ds)):
        c=cost(gs,ds)
        print(f"  *** cost {c}: break {gs}, drop {ds}  ({time.time()-t0:.0f}s)", flush=True)
        if best is None or c<best[0]: best=(c,gs,ds)
print("BEST:",best,f"({time.time()-t0:.0f}s)")
if best: json.dump({'cost':best[0],'gates':list(best[1]),'drop':list(best[2])},
                   open(os.path.join(HERE,'data','gmp33.json'),'w'))
