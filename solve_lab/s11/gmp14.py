"""The MAXIMAL continuous move set: all 1,726 live free inputs at once.

If this system is consistent over GF(p), the instance has a state with every check zero mod p
-- and then the p-quantised handles absorb the residues exactly.  If it is inconsistent, no
continuous move from this base can do it and only bit flips remain, which is itself a definitive
statement (every earlier obstruction proof was about a far smaller move set).
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
F=failing(bd)
rows=set(F)
for u,d in cols.items(): rows |= set(d)
rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
knobs=sorted(cols)
print(f"system {len(rowl)} rows x {len(knobs)} knobs; failing {len(F)} {F}")
nz=sum(len(d) for d in cols.values()); print(f"nonzeros {nz}")
# sparse column-major elimination with [I] tracking for certificates
R=[dict() for _ in rowl]
for j,u in enumerate(knobs):
    for a,x in cols[u].items(): R[idx[a]][j]=x%P
b=[(-bd[a])%P for a in rowl]
E=[{i:1} for i in range(len(rowl))]        # left-multiplier bookkeeping
m=len(rowl); n=len(knobs)
used=[False]*m; piv={}
t0=time.time()
for c in range(n):
    cand=[i for i in range(m) if not used[i] and c in R[i]]
    if not cand: continue
    i=min(cand, key=lambda i: len(R[i]))
    used[i]=True; piv[c]=i
    inv=pow(R[i][c],-1,P)
    R[i]={k:x*inv%P for k,x in R[i].items()}
    E[i]={k:x*inv%P for k,x in E[i].items()}
    b[i]=b[i]*inv%P
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
    if c%200==0: print(f"   col {c}/{n} rank {len(piv)} ({time.time()-t0:.0f}s)", flush=True)
print(f"rank {len(piv)} of min({m},{n})  ({time.time()-t0:.0f}s)")
bad=[i for i in range(m) if not R[i] and b[i]]
print("inconsistent rows:",len(bad))
if bad:
    for i in bad[:6]:
        sup=[(rowl[k],v) for k,v in E[i].items()]
        print(f"   certificate support {len(sup)}: {[rowl[k] for k in E[i]]}")
else:
    x=[0]*n
    for c,i in piv.items(): x[c]=b[i]
    v=list(base)
    for j,u in enumerate(knobs): v[u]=(v[u]+x[j])%P
    forwardp(v)
    F2=failing(allchk(v))
    print(f"APPLIED -> failing checks mod p = {len(F2)} {F2[:20]}")
    json.dump([int(t) for t in v], open(os.path.join(HERE,'data','gmp14_state.json'),'w'))
    print("saved data/gmp14_state.json")
