"""Break a36244 (4 equations) to buy the missing knob, then solve the mod-p system.

The continuous deficit is exactly 1, and a36244 is a gate atom lying in only FOUR equations whose
freed output x3432 moves a25676, a42245, a29539 and a40826 -- the deficit rows.  Freeze x3432
(so the gate stays broken and x3432 becomes a genuine knob), rebuild the exact response matrix
with that freeze in place, and solve.  If it closes, the only failing equations are a36244's four.
"""
import sys, os, json, time, pickle, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
from gmp26 import forwardp_frozen
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=[u for u in range(L.NVARS) if u not in L.definer]
GATE_BREAK=int(sys.argv[1]) if len(sys.argv)>1 else 36244
T=L.atom_out[GATE_BREAK][1]
FROZEN={T}
KNOBS=sorted(set(FREE)|{T})
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp_frozen(base,FROZEN)
bd={a:evalp(L.polys[a],base) for a in CHK}
F=[a for a in CHK if bd[a]]
print(f"breaking a{GATE_BREAK} (cost {len(L.atom2eq.get(GATE_BREAK,{}))} eqs), knob x{T}")
print(f"  failing checks mod p: {len(F)} {F}")
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
t0=time.time(); cols={}
for i,u in enumerate(KNOBS):
    if u in BITS: continue
    v=list(base); v[u]=(v[u]+1)%P; forwardp_frozen(v,FROZEN)
    d={}
    for a in CHK:
        x=evalp(L.polys[a],v)
        if x!=bd[a]: d[a]=(x-bd[a])%P
    if d: cols[u]=d
    if i%800==0: print(f"   {i}/{len(KNOBS)} live={len(cols)} ({time.time()-t0:.0f}s)", flush=True)
print(f"live non-bit knobs: {len(cols)}  ({time.time()-t0:.0f}s)")
rows=set(F)
for u in cols: rows |= set(cols[u])
rowl=sorted(rows); idx={a:i for i,a in enumerate(rowl)}
knobs=sorted(cols); n=len(knobs); m=len(rowl)
print(f"system {m} x {n}")
R=[dict() for _ in rowl]
for j,u in enumerate(knobs):
    for a,x in cols[u].items(): R[idx[a]][j]=x%P
b=[(-bd[a])%P for a in rowl]
used=[False]*m; piv={}
for c in range(n):
    cand=[i for i in range(m) if not used[i] and c in R[i]]
    if not cand: continue
    i=min(cand,key=lambda i:len(R[i])); used[i]=True; piv[c]=i
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
bad=[i for i in range(m) if not R[i] and b[i]]
print(f"rank {len(piv)}; inconsistent rows {len(bad)} -> {[rowl[i] for i in bad][:10]}")
if not bad:
    x=[0]*n
    for c,i in piv.items(): x[c]=b[i]
    v=list(base)
    for j,u in enumerate(knobs): v[u]=(v[u]+x[j])%P
    forwardp_frozen(v,FROZEN)
    F2=[a for a in CHK if evalp(L.polys[a],v)]
    AV=[evalp(L.polys[a],v) for a in range(L.NA)]
    E=[e for e in range(L.NEQ) if sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())%P]
    print(f"APPLIED -> failing checks mod p = {len(F2)} {F2[:10]}; equations nonzero mod p = {len(E)}")
    json.dump([int(t) for t in v], open(os.path.join(HERE,'data','gmp27_state.json'),'w'))
    print("saved data/gmp27_state.json")
