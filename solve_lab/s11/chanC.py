"""Is the UNEXPLORED third channel (U=1, V=0, x34606) mod-p solvable?"""
import sys, os, json, time, pickle, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
D=pickle.load(open(sys.argv[1],'rb'))
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
rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
m=len(rowl); n=len(knobs)
print(f"channel system {m} rows x {n} non-bit knobs; failing checks {len(F)} {F}")
R=[dict() for _ in rowl]
for j,u in enumerate(knobs):
    for a,x in cols[u].items(): R[idx[a]][j]=x%P
b=[(-bd[a])%P for a in rowl]
E=[{i:1} for i in range(m)]
used=[False]*m; piv={}
for c in range(n):
    cand=[i for i in range(m) if not used[i] and c in R[i]]
    if not cand: continue
    i=min(cand,key=lambda i:len(R[i])); used[i]=True; piv[c]=i
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
bad=[i for i in range(m) if not R[i] and b[i]]
print(f"rank {len(piv)}; INCONSISTENT rows {len(bad)}")
seen=set()
for i in bad:
    key=tuple(sorted(rowl[k] for k in E[i]))
    if key in seen: continue
    seen.add(key)
    print(f"   certificate ({len(key)} rows): {list(key)}   failing among them: {[a for a in key if bd.get(a)]}")
if not bad:
    print("*** CHANNEL IS MOD-P SOLVABLE ***")
    x=[0]*n
    for c,i in piv.items(): x[c]=b[i]
    v=list(base)
    for j,u in enumerate(knobs): v[u]=(v[u]+x[j])%P
    forwardp(v)
    F2=[a for a in CHK if evalp(L.polys[a],v)]
    print("applied -> failing checks:",len(F2),F2[:10])
    json.dump([int(t) for t in v], open(os.path.join(HERE,'data','chanC_solved.json'),'w'))
