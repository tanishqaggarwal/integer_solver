#!/usr/bin/env python3
"""MAIN quadratic-core solver (agentC).
Controls a=x_12186, b=x_14853, c=x_16742 (freeze x_24908 cone so x_3558=x_24908-c moves only via c).
Build compensator-closure handle set + touched-sat constraints; null space; project onto (a,b,c);
solve S(t)=T(t)=0 on the achievable subspace using the QUADRATIC structure (Tonelli-Shanks)."""
import json, time, sys
from agentC_common import (p, gates, order, definer, gcode, forward, partial_forward, downstream_ks,
                           val, freeinp, ns, lines, eqcode, eqvars, load_best, CORE, posof,
                           NVARS, pinned, rootcode_of, inv, is_qr, sqrt_mod, C1, C2)
from collections import defaultdict

gate_defs = {t: (rhs, vids) for t, rhs, vids in gates}
def freecone(root):
    seen=set(); leaves=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x][1]: st.append(u)
        elif x in freeinp: leaves.add(x)
    return leaves

best=load_best(); forward()
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
# deep cones (to exclude from compensators so they don't move deep vars)
DEEP=[3558,29322,33469,27713,1326]
deepcone=set()
for d in DEEP: deepcone|=freecone(d)
deepcone|= (freecone(35389)|freecone(6671))
a_id,b_id,c_id=12186,14853,16742
CONTROLS=[a_id,b_id,c_id]
# --- compensator closure ---
# handles start with controls; constraints = sat-eqs touched; add free inputs of broken eqs as
# compensators (excluding deepcone); iterate to fixpoint but bounded.
def affected_eqs(h):
    ks=downstream_ks(h)
    aff=set(eqbyvar.get(h,()))
    for k in ks: aff|=eqbyvar.get(order[k],set())
    return aff
H=set(CONTROLS)
cons=set()   # constraint equations (currently-satisfied)
for _hop in range(6):
    newcons=set()
    for h in H: newcons|=(affected_eqs(h)-F0)
    cons=newcons
    # compensators = free inputs in constraint eqs, minus deepcone, minus controls
    comp=set()
    for i in cons: comp|=(eqvars[i]&freeinp)
    comp-=deepcone
    comp|=set(CONTROLS)
    if comp<=H:
        break
    H=comp
    print(f"  hop {_hop}: handles={len(H)}, constraints={len(cons)}", flush=True)
H=sorted(H)
Hidx={h:i for i,h in enumerate(H)}
NH=len(H)
cons=sorted(cons)
print(f"FINAL handle set: {NH}; constraints (touched-sat): {len(cons)}", flush=True)
# sanity: controls in H, deepcone∩H == controls only
assert all(cc in Hidx for cc in CONTROLS)
overlap=(set(H)&deepcone)-set(CONTROLS)
print(f"handle∩deepcone beyond controls: {len(overlap)} (should be 0): {sorted(overlap)[:10]}")

# --- Build Jacobian J (constraints x handles) mod p via finite diff ---
base_root={i:eval(rootcode_of(i),ns)%p for i in cons}
nz=[i for i in cons if base_root[i]!=0]
print(f"constraints with nonzero root (unexpected): {len(nz)}")
dks={h:downstream_ks(h) for h in H}
t0=time.time()
Jcols=[]
consset=set(cons)
aff_h={}
for h in H:
    a=set(eqbyvar.get(h,()))&consset
    for k in dks[h]: a|=(eqbyvar.get(order[k],())&consset)
    aff_h[h]=sorted(a)
for hi,h in enumerate(H):
    o=val[h]; val[h]=o+1; partial_forward(dks[h])
    col={}
    for i in aff_h[h]:
        d=(eval(rootcode_of(i),ns)-base_root[i])%p
        if d: col[i]=d
    Jcols.append(col); val[h]=o; partial_forward(dks[h])
    if hi%100==0: print(f"  jac {hi}/{NH} ({time.time()-t0:.0f}s)", flush=True)
print(f"jacobian built {time.time()-t0:.0f}s", flush=True)
rows_used=sorted(set().union(*[set(c) for c in Jcols])) if Jcols else []
ridx={r:i for i,r in enumerate(rows_used)}
M=len(rows_used)
Jm=[[0]*NH for _ in range(M)]
for hj,col in enumerate(Jcols):
    for r,vv in col.items(): Jm[ridx[r]][hj]=vv%p
print(f"active constraint rows: {M}", flush=True)

def rref(Mrows, ncol):
    Mm=[row[:] for row in Mrows]; m=len(Mm); r=0; piv=[]
    for c in range(ncol):
        sel=None
        for i in range(r,m):
            if Mm[i][c]%p: sel=i;break
        if sel is None: continue
        Mm[r],Mm[sel]=Mm[sel],Mm[r]; iv=inv(Mm[r][c]); Mm[r]=[(x*iv)%p for x in Mm[r]]
        for i in range(m):
            if i!=r and Mm[i][c]%p:
                f=Mm[i][c]; Mm[i]=[(Mm[i][k]-f*Mm[r][k])%p for k in range(ncol)]
        piv.append(c); r+=1
        if r>=m: break
    return Mm[:r], piv

Rrows,piv=rref(Jm,NH)
rank=len(piv); nullity=NH-rank
print(f"J rank={rank}, nullity={nullity} over {NH} handles", flush=True)
pivset=set(piv); free_cols=[c for c in range(NH) if c not in pivset]
Null=[]
for f in free_cols:
    n=[0]*NH; n[f]=1
    for ri,c in enumerate(piv): n[c]=(-Rrows[ri][f])%p
    Null.append(n)
ia,ib,ic=Hidx[a_id],Hidx[b_id],Hidx[c_id]
Proj=[[n[ia]%p,n[ib]%p,n[ic]%p] for n in Null]
# rank of achievable (da,db,dc)
Pr,Pp=rref(Proj,3)
dimV=len(Pp)
print(f"achievable (da,db,dc) subspace dim = {dimV}", flush=True)
# select dimV null vectors with independent projections
basis_idx=[]; chosen=[]
Ptmp=[]
for j,pj in enumerate(Proj):
    trial=Ptmp+[pj]
    if len(rref(trial,3)[1])>len(rref(Ptmp,3)[1]):
        Ptmp.append(pj); basis_idx.append(j)
        if len(basis_idx)==dimV: break
projbasis=[Proj[j] for j in basis_idx]     # dimV vectors in F_p^3
nullbasis=[Null[j] for j in basis_idx]     # corresponding full moves
print(f"projection basis: {projbasis}", flush=True)

# --- deep var residues as affine functions of t (params of achievable subspace) ---
# current residues
xa=val[a_id]%p; xb=val[b_id]%p; xc=val[c_id]%p; e=val[24908]%p
K33=(C1 + 97553848499418123410591666447050222001188385549510401465815187079080512838891)%p
x29322=val[29322]%p; x1326=val[1326]%p; x27713=val[27713]%p; x33469=val[33469]%p; x3558=val[3558]%p
# verify constant K33 and identities
assert (xb-xa)%p==x29322
assert (xa-C1)%p==x1326
assert (C2+xc)%p==x27713
assert (xa+xb+K33)%p==x33469, ((xa+xb+K33)%p, x33469)
assert (e-xc)%p==x3558
# (da,db,dc)(t) = sum_j t_j * projbasis[j]
def deep_of_t(tv):
    da=sum(tv[j]*projbasis[j][0] for j in range(dimV))%p
    db=sum(tv[j]*projbasis[j][1] for j in range(dimV))%p
    dc=sum(tv[j]*projbasis[j][2] for j in range(dimV))%p
    A=(xa+da)%p; B=(xb+db)%p; Cc=(xc+dc)%p
    X29322=(B-A)%p; X1326=(A-C1)%p; X27713=(C2+Cc)%p; X33469=(A+B+K33)%p; X3558=(e-Cc)%p
    S=(X33469*X29322*X29322 - X3558*X3558)%p
    T=(X27713*X29322 - X3558*X1326)%p
    return S,T
# sanity at t=0
S0,T0=deep_of_t([0]*dimV)
print(f"S(0)={S0} (should={val[35389]%p}); T(0)={T0} (should={val[6671]%p})")

json.dump({'H':H,'dimV':dimV,'projbasis':projbasis,'nullbasis':nullbasis,
           'xa':xa,'xb':xb,'xc':xc,'e':e,'K33':K33},
          open('agentC_Vdata.json','w'))
print("saved agentC_Vdata.json")
