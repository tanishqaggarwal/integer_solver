"""The obstruction with the bits removed.

Continuously the bits are frozen: each bit is the ONLY knob touching its own pin row, so any
linear system containing the pin rows forces the bit coefficient to zero.  Dropping the bits (and
their pin rows) leaves the genuine continuous system, whose certificates should be small and
readable rather than 283 rows of bookkeeping.
"""
import sys, os, json, time, pickle, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
from gmp7 import allchk, failing
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
print(f"non-bit system: {len(rowl)} rows x {len(knobs)} knobs; failing {F}")
m=len(rowl); n=len(knobs)
R=[dict() for _ in rowl]
for j,u in enumerate(knobs):
    for a,x in cols[u].items(): R[idx[a]][j]=x%P
b=[(-bd[a])%P for a in rowl]
E=[{i:1} for i in range(m)]
used=[False]*m; piv={}
for c in range(n):
    cand=[i for i in range(m) if not used[i] and c in R[i]]
    if not cand: continue
    i=min(cand, key=lambda i: len(R[i]))
    used[i]=True; piv[c]=i
    inv=pow(R[i][c],-1,P)
    R[i]={k:x*inv%P for k,x in R[i].items()}; E[i]={k:x*inv%P for k,x in E[i].items()}; b[i]=b[i]*inv%P
    Ri=R[i]; Ei=E[i]; bi=b[i]
    for k in cand:
        if k==i: continue
        f=R[k].get(c)
        if not f: continue
        Rk=R[k]
        for kk,x in Ri.items():
            nv=(Rk.get(kk,0)-f*x)%P
            if nv: Rk[kk]=nv
            elif kk in Rk: del Rk[kk]
        Ek=E[k]
        for kk,x in Ei.items():
            nv=(Ek.get(kk,0)-f*x)%P
            if nv: Ek[kk]=nv
            elif kk in Ek: del Ek[kk]
        b[k]=(b[k]-f*bi)%P
print("rank",len(piv))
bad=[i for i in range(m) if not R[i] and b[i]]
print("inconsistent rows:",len(bad))
seen=set()
for i in bad:
    sup=sorted(E[i])
    key=tuple(rowl[k] for k in sup)
    if key in seen: continue
    seen.add(key)
    print(f"  certificate ({len(sup)} rows): {[rowl[k] for k in sup]}")
    print(f"     failing among them: {[rowl[k] for k in sup if bd[rowl[k]]]}")
