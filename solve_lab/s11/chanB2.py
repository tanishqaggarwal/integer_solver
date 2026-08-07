"""Channel B's own certificates -> invariants -> which bits move them.

Channel B is the closest to solvable: 3 failing checks, only 4 obstruction directions, and
a26723 appears in none of them (so it is fixable).  Derive its certificates properly and repeat
the dependency analysis that made channel A's inv5 enumerable.
"""
import sys, os, json, time, pickle, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
D=pickle.load(open(os.path.join(HERE,'data','resp_B.pkl'),'rb'))
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
R=[dict() for _ in rowl]
for j,u in enumerate(knobs):
    for a,x in cols[u].items(): R[idx[a]][j]=x%P
b=[(-bd[a])%P for a in rowl]
E=[{i:1} for i in range(m)]
used=[False]*m
for c in range(n):
    cand=[i for i in range(m) if not used[i] and c in R[i]]
    if not cand: continue
    i=min(cand,key=lambda i:len(R[i])); used[i]=True
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
CERT=[{rowl[k]:v for k,v in E[i].items()} for i in bad]
print(f"channel B: {len(CERT)} certificates, supports {[len(c) for c in CERT]}")
json.dump([{str(k):str(v) for k,v in c.items()} for c in CERT],
          open(os.path.join(HERE,'data','certsB.json'),'w'))
# invariants and their bit dependency, from the 39,018 message {x91, x490}
from bits5 import msg, real, tree
def inv(v): return tuple(sum(y*evalp(L.polys[a],v) for a,y in c.items())%P for c in CERT)
ANCH={91,490}
v0=msg(ANCH); I0=inv(v0)
print("invariants at the 39,018 message:", [('ZERO' if not x else str(x)[:14]+'..') for x in I0])
dep=collections.defaultdict(list)
t0=time.time()
for bb in real:
    S=(ANCH-{bb}) if bb in ANCH else (ANCH|{bb})
    if not S: continue
    I=inv(msg(S))
    for j in range(len(CERT)):
        if I[j]!=I0[j]: dep[j].append(bb)
for j in range(len(CERT)):
    d=dep[j]; byq=collections.Counter(tree.get(x,'?') for x in d)
    print(f"   invB{j}: moved by {len(d):3d} bits {dict(byq)}")
    if len(d)<=40: print(f"       {sorted(d)}")
print(f"({time.time()-t0:.0f}s)")
json.dump({str(j):sorted(dep[j]) for j in dep}, open(os.path.join(HERE,'data','chanB2.json'),'w'))
