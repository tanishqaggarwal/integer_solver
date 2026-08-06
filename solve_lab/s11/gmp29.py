"""Cheapest SET of gates to break so the mod-p system closes.

Breaking one gate is not enough: the residual has to land inside the enlarged column space, not
merely break some certificate.  So measure the freed-output column of every cheap gate, add them
all, check consistency, then shrink to a minimal-cost subset.  Cost = |equations containing the
broken gate atoms|, which is exactly the score penalty.
"""
import sys, os, json, time, pickle, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
from gmp26 import forwardp_frozen
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
D=pickle.load(open(os.path.join(HERE,'data','resp_modp.pkl'),'rb'))
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
bd={a:evalp(L.polys[a],base) for a in CHK}
F=[a for a in CHK if bd[a]]
knobs=sorted(u for u in D['cols'] if u not in BITS)
cache=os.path.join(HERE,'data','gmp29_cols.pkl')
if os.path.exists(cache):
    Z=pickle.load(open(cache,'rb')); C=Z['C']; G=Z['G']
else:
    t0=time.time(); C={}
    for u in knobs:
        v=list(base); v[u]=(v[u]+1)%P; forwardp_frozen(v,set())
        d={a:(evalp(L.polys[a],v)-bd[a])%P for a in CHK}
        d={a:x for a,x in d.items() if x}
        if d: C[u]=d
    print(f"knob responses: {len(C)} ({time.time()-t0:.0f}s)", flush=True)
    GATES=[a for a in range(L.NA) if L.atom_out.get(a) is not None]
    cost={a:len(L.atom2eq.get(a,{})) for a in GATES}
    cands=[g for g in sorted(GATES,key=lambda a:cost[a]) if cost[g]<=8]
    t0=time.time(); G={}
    for g in cands:
        t=L.atom_out[g][1]
        v=list(base); v[t]=(v[t]+1)%P; forwardp_frozen(v,{t})
        d={a:(evalp(L.polys[a],v)-bd[a])%P for a in CHK}
        d={a:x for a,x in d.items() if x}
        if d: G[g]=d
    print(f"gate columns: {len(G)} ({time.time()-t0:.0f}s)", flush=True)
    pickle.dump({'C':C,'G':G}, open(cache,'wb'))
COST={g:len(L.atom2eq.get(g,{})) for g in G}
def consistent(extra):
    cols=[C[u] for u in knobs if u in C]+[G[g] for g in extra]
    rows=set(F)
    for d in cols: rows |= set(d)
    rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
    m=len(rowl); n=len(cols)
    R=[dict() for _ in rowl]
    for j,d in enumerate(cols):
        for a,x in d.items(): R[idx[a]][j]=x
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
allg=sorted(G, key=lambda g:COST[g])
print(f"candidates {len(allg)}; failing checks {F}")
t0=time.time()
print("consistent with NO gates broken:", consistent([]), flush=True)
print("consistent with ALL candidate gates broken:", consistent(allg), f"({time.time()-t0:.0f}s)", flush=True)
