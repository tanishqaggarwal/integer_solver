"""Cheap CLUSTERS, tested against the mod-p solver.

A cluster is a set of atoms whose equations overlap so heavily that breaking all of them costs
only |union of their equations|.  The 39,026 checkpoint is exactly such a cluster, of cost 7.
Enumerate every cluster of cost <= 6 (all atoms whose equation-set is contained in the union),
freeze the gate outputs it frees, allow its checks to be nonzero, and test whether everything
else can be made zero mod p.
"""
import sys, os, json, time, pickle, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
from gmp26 import forwardp_frozen
P=L.P; sys.set_int_max_str_digits(400000)
CHKS=set(a for a in range(L.NA) if L.atom_out.get(a) is None)
EQ={a:frozenset(L.atom2eq.get(a,{})) for a in range(L.NA)}
cheap=[a for a in range(L.NA) if 0<len(EQ[a])<=6]
print(f"atoms in <=6 equations: {len(cheap)}  ({sum(1 for a in cheap if a not in CHKS)} gates)")
byeq=collections.defaultdict(list)
for a in cheap:
    for e in EQ[a]: byeq[e].append(a)
def cluster_of(u):
    cand=set(byeq[next(iter(u))])
    for e in u: cand |= set(byeq[e])
    return tuple(sorted(a for a in cand if EQ[a]<=u))
U={}
for a in cheap:
    u=EQ[a]
    if len(u)<=6: U[u]=None
for a,b in itertools.combinations(cheap,2):
    u=EQ[a]|EQ[b]
    if len(u)<=6: U[u]=None
print(f"candidate equation-unions of size <=6: {len(U)}")
CL={}
for u in U:
    c=cluster_of(u)
    if c and (c not in CL or len(u)<CL[c]): CL[c]=len(u)
print(f"distinct clusters: {len(CL)}")
top=sorted(CL.items(), key=lambda z:(z[1], -sum(1 for a in z[0] if a not in CHKS)))
for c,cost in top[:12]:
    print(f"   cost {cost}: {len(c)} atoms ({sum(1 for a in c if a not in CHKS)} gates) {c[:8]}")
# ---- test them
Z=pickle.load(open(os.path.join(HERE,'data','gmp29_cols.pkl'),'rb')); C=Z['C']
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
bd={a:evalp(L.polys[a],base) for a in CHKS}
F=[a for a in CHKS if bd[a]]
knobcols=[C[u] for u in sorted(C)]
def colfor(t, frozen):
    v=list(base); v[t]=(v[t]+1)%P; forwardp_frozen(v,frozen)
    d={a:(evalp(L.polys[a],v)-bd[a])%P for a in CHKS}
    return {a:x for a,x in d.items() if x}
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
t0=time.time(); hits=[]
for ci,(c,cost) in enumerate(top):
    gates=[a for a in c if a not in CHKS]
    drops=[a for a in c if a in CHKS]
    frozen={L.atom_out[a][1] for a in gates}
    cols=knobcols+[colfor(t,frozen) for t in frozen]
    if consistent(cols, set(drops)):
        hits.append((cost,c)); print(f"  *** cost {cost} cluster CLOSES: {c}", flush=True)
    if ci%150==0: print(f"   {ci}/{len(top)} ({time.time()-t0:.0f}s)", flush=True)
print(f"clusters that close the system: {len(hits)} ({time.time()-t0:.0f}s)")
for cost,c in sorted(hits)[:10]: print("   cost",cost,c)
