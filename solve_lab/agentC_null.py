#!/usr/bin/env python3
"""Compute null space of the touched-satisfied equations (roots mod p) over the 97 deep-cone handles.
Measure: (1) rank/nullity, (2) whether x_29322, x_3558 residues are movable in the null space
(regime-1 linear feasibility), (3) achievable (a,b,c,e) subspace dimension for the quadratic solve."""
import json, time, sys
from agentC_common import (p, gates, order, definer, gcode, forward, partial_forward, downstream_ks,
                           val, freeinp, anc, ns, lines, eqcode, eqvars, load_best, CORE, posof,
                           NVARS, pinned, rootcode_of, inv, is_qr, sqrt_mod)
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

best = load_best(); forward()
resp = json.load(open('agentC_resp.json'))
H = resp['handles']              # 97 handles
Hidx = {h:i for i,h in enumerate(H)}
NH = len(H)
Hset = set(H)
# baseline failing (core)
F0 = set(i for i in range(len(lines)) if eval(eqcode[i], ns)!=0)
# touched-satisfied equations
TS = sorted(i for i in range(len(lines)) if (eqvars[i]&Hset) and i not in F0)
print(f"handles={NH}, touched-satisfied eqs={len(TS)}", flush=True)
# eqs affected by each handle (via downstream gates too) for fast Jacobian
eqbyvar = defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
# precompute downstream ks + affected TS eqs per handle
dks={}; aff={}
TSset=set(TS)
for h in H:
    ks = downstream_ks(h)
    dks[h]=ks
    a=set(eqbyvar.get(h,())) & TSset
    for k in ks: a |= (eqbyvar.get(order[k],()) & TSset)
    aff[h]=sorted(a)
# baseline roots mod p for TS
base_root = {i: eval(rootcode_of(i), ns)%p for i in TS}
nz = [i for i in TS if base_root[i]!=0]
print(f"TS eqs with nonzero root mod p (unexpected): {len(nz)} (first {nz[:5]})", flush=True)
# Build Jacobian J: rows=TS, cols=H, entry = d(root_i)/d(h) mod p
# store sparse: col-major dict
t0=time.time()
Jcols=[]  # list over handles of dict{row_i: val}
for hi,h in enumerate(H):
    o=val[h]; val[h]=o+1
    partial_forward(dks[h])
    col={}
    for i in aff[h]:
        d=(eval(rootcode_of(i), ns)-base_root[i])%p
        if d: col[i]=d
    Jcols.append(col)
    val[h]=o; partial_forward(dks[h])
    if hi%25==0: print(f"  jac {hi}/{NH} ({time.time()-t0:.0f}s)", flush=True)
print(f"jacobian built in {time.time()-t0:.0f}s", flush=True)

# rows that appear
rows_used = sorted(set().union(*[set(c) for c in Jcols])) if Jcols else []
ridx={r:i for i,r in enumerate(rows_used)}
M=len(rows_used)
print(f"active constraint rows: {M}", flush=True)
# Dense matrix Jm  (M x NH)
Jm=[[0]*NH for _ in range(M)]
for hj,col in enumerate(Jcols):
    for r,vv in col.items():
        Jm[ridx[r]][hj]=vv%p

# --- deep var linear responses to handles (exact, for a,b,c,e and the 5 deep residues) ---
# a=x_12186 idx, b=x_14853, c=x_16742 are handles. e=x_24908 via cone response.
# measure response of x_24908 and the 5 deep vars to each handle (mod p) reuse from resp? resp had deep vars.
# recompute x_24908 response:
base_e = val[24908]%p
e_resp={}
for h in H:
    o=val[h]; val[h]=o+1; partial_forward(dks[h])
    d=(val[24908]%p-base_e)%p
    if d: e_resp[h]=d
    val[h]=o; partial_forward(dks[h])
# deep residue responses (from resp json, keyed by deep var id)
DEEP=[3558,29322,33469,27713,1326]
Dresp={d:{} for d in DEEP}
for h in H:
    hr=resp['resp'].get(str(h),{})
    for d in DEEP:
        if str(d) in hr: Dresp[d][Hidx[h]]=int(hr[str(d)])%p

def rref(M, ncol):
    """Return (rref matrix rows, pivot cols). M is list of rows length ncol."""
    M=[row[:] for row in M]; m=len(M); r=0; piv=[]
    for c in range(ncol):
        sel=None
        for i in range(r,m):
            if M[i][c]%p: sel=i;break
        if sel is None: continue
        M[r],M[sel]=M[sel],M[r]; iv=inv(M[r][c]); M[r]=[(x*iv)%p for x in M[r]]
        for i in range(m):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(M[i][k]-f*M[r][k])%p for k in range(ncol)]
        piv.append(c); r+=1
        if r>=m: break
    return M[:r], piv

Rrows,piv = rref(Jm, NH)
rank=len(piv); nullity=NH-rank
print(f"J rank={rank}, nullity={nullity} over {NH} handles", flush=True)
# null space basis: free cols = non-pivot
pivset=set(piv)
free_cols=[c for c in range(NH) if c not in pivset]
# for each free col f, basis vector n with n[f]=1, n[piv]= -Rrows[row][f]
Null=[]
for f in free_cols:
    n=[0]*NH; n[f]=1
    for ri,c in enumerate(piv):
        n[c]=(-Rrows[ri][f])%p
    Null.append(n)
print(f"null basis vectors: {len(Null)}", flush=True)

# achievable deep residue moves: for each null vector, (d3558,d29322,d33469,d27713,d1326)
def apply_resp(Dr, n):
    return sum(Dr.get(j,0)*n[j] for j in range(NH))%p
Amat=[]  # rows = null vectors, cols = 5 deep deltas
for n in Null:
    Amat.append([apply_resp(Dresp[d], n) for d in DEEP])
# dimension of achievable deep-residue space
Ar,Ap = rref(Amat, 5)
print(f"achievable deep-residue subspace dim = {len(Ap)} (of 5): pivot deep vars {[DEEP[c] for c in Ap]}", flush=True)

# Regime-1 linear feasibility: exists n in null with d(x_29322) = -x_29322 and d(x_3558) = -x_3558 ?
# d(x_29322)=apply_resp(Dresp[29322],n); want == (-val[29322])%p ; similarly x_3558.
# Solve over null coefficients: sum_j alpha_j * (Amat[j] restricted to [29322,3558]) = targets
i29=DEEP.index(29322); i35=DEEP.index(3558)
tgt=[(-val[29322])%p,(-val[3558])%p]
# system: for each null vec j, coeff alpha_j; rows = 2 (the two deep vars)
Csys=[[Amat[j][i29] for j in range(len(Null))],
      [Amat[j][i35] for j in range(len(Null))]]
# solve Csys * alpha = tgt (2 equations, len(Null) unknowns) -> augmented rref
aug=[Csys[0][:]+[tgt[0]], Csys[1][:]+[tgt[1]]]
ar,ap=rref(aug, len(Null)+1)
# consistent if no pivot in last col
regime1_ok = (len(Null) not in ap)
print(f"regime-1 (x_29322->0 AND x_3558->0) linear-feasible in null space: {regime1_ok}", flush=True)

json.dump({'H':H,'rank':rank,'nullity':nullity,
           'achievable_dim':len(Ap),
           'regime1_ok':regime1_ok}, open('agentC_null_summary.json','w'))
print("saved agentC_null_summary.json")
