"""Which gate, if broken, supplies the missing knob direction?

The obstruction is a left vector y with y^T J = 0 over every continuous knob and y . rhs != 0.
Breaking gate atom g frees its output t; the new knob helps iff  y . J_t != 0.  That is a single
dot product per candidate, so every gate can be tested -- and the cost of the winner is just the
number of equations containing g.
"""
import sys, os, json, time, pickle, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
from gmp26 import forwardp_frozen
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
D=pickle.load(open(os.path.join(HERE,'data','resp_modp.pkl'),'rb'))
base0=D['base']; cols=D['cols']
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
# rebuild the base with the two free checks cleared (matches gmp16_base)
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
bd={a:evalp(L.polys[a],base) for a in CHK}
F=[a for a in CHK if bd[a]]
print("failing checks:",F)
knobs=sorted(u for u in cols if u not in BITS)
# exact responses at THIS base
t0=time.time(); C={}
for u in knobs:
    v=list(base); v[u]=(v[u]+1)%P; forwardp_frozen(v,set())
    d={a:(evalp(L.polys[a],v)-bd[a])%P for a in CHK}
    d={a:x for a,x in d.items() if x}
    if d: C[u]=d
print(f"re-measured {len(C)} knob responses at the base ({time.time()-t0:.0f}s)", flush=True)
rows=set(F)
for u in C: rows |= set(C[u])
rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
kl=sorted(C); n=len(kl); m=len(rowl)
R=[dict() for _ in rowl]
for j,u in enumerate(kl):
    for a,x in C[u].items(): R[idx[a]][j]=x
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
print(f"system {m}x{n}; inconsistent rows {len(bad)}")
Y=[{rowl[k]:v for k,v in E[i].items()} for i in bad]
print("certificate supports:", [len(y) for y in Y])
json.dump([{str(k):str(v) for k,v in y.items()} for y in Y], open(os.path.join(HERE,'data','certs.json'),'w'))
GATES=[a for a in range(L.NA) if L.atom_out.get(a) is not None]
cost={a:len(L.atom2eq.get(a,{})) for a in GATES}
cands=sorted(GATES,key=lambda a:cost[a])
cands=[g for g in cands if cost[g]<=8]
print(f"testing {len(cands)} gate atoms with cost <= 8", flush=True)
t0=time.time(); hits=[]
for i,g in enumerate(cands):
    t=L.atom_out[g][1]
    v=list(base); v[t]=(v[t]+1)%P; forwardp_frozen(v,{t})
    col={a:(evalp(L.polys[a],v)-bd[a])%P for a in CHK}
    col={a:x for a,x in col.items() if x}
    good=[j for j,y in enumerate(Y) if sum(y.get(a,0)*x for a,x in col.items())%P]
    if good: hits.append((cost[g],g,t,good,len(col)))
    if i%200==0: print(f"   {i}/{len(cands)} hits={len(hits)} ({time.time()-t0:.0f}s)", flush=True)
hits.sort()
print(f"gates that break at least one certificate: {len(hits)}")
for c,g,t,good,nc in hits[:25]:
    print(f"   a{g} cost={c} eqs frees x{t}: breaks certificates {good}  moves {nc} checks")
json.dump([[c,g,t,good,nc] for c,g,t,good,nc in hits], open(os.path.join(HERE,'data','gmp28.json'),'w'))
