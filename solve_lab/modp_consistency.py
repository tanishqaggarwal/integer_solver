#!/usr/bin/env python3
import heal_harness as H
import time
from collections import defaultdict
p=H.p
def inv(a): return pow(a%p,p-2,p)
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
ns={'v':H.val,'__builtins__':{}}
F16=set([697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431])
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
eq_free=[]
for i in range(len(H.eqcode)):
    s=set()
    for w in H.eqvars[i]:
        if w in H.freeinp: s.add(w)
        s|=H.anc.get(w,set())
    eq_free.append(s&H.freeinp)
free_eqs=defaultdict(set)
for i in range(len(H.eqcode)):
    for w in eq_free[i]: free_eqs[w].add(i)
pool=set().union(*[eq_free[i] for i in F16])
eqs=set(F16)
CAP=900
for _ in range(50):
    ne=set()
    for w in pool: ne|=free_eqs[w]
    nf=set()
    for i in ne: nf|=eq_free[i]
    if len(pool|nf)>CAP:
        for w in sorted(nf-pool):
            if len(pool)>=CAP: break
            pool.add(w)
        eqs|=ne; break
    if ne<=eqs and nf<=pool: break
    eqs|=pool and ne; eqs|=ne; pool|=nf
poolL=sorted(pool); Et=sorted(eqs)
print(f"pool={len(poolL)} frees, touched eqs={len(Et)}",flush=True)
t0=time.time()
r0={i:eval(H.eqcode[i],ns)%p for i in Et}
base=H.val[:]
col={w:j for j,w in enumerate(poolL)}; nc=len(poolL)
# rows: eq -> {col: coef}
rowdict=defaultdict(dict)
for wi,w in enumerate(poolL):
    H.val[w]+=1
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
    for i in Et:
        nv=eval(H.eqcode[i],ns)%p
        d=(nv-r0[i])%p
        if d: rowdict[i][col[w]]=d
    H.val[w]=base[w]
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
    if wi%150==0: print(f"  jac {wi}/{nc} ({time.time()-t0:.0f}s)",flush=True)
# Build augmented system: for F16 rhs=-r0 (want 0), for others rhs=0 (keep 0)
A=[]; ridx=[]
for i in Et:
    row=[0]*(nc+1)
    for j,c in rowdict[i].items(): row[j]=c
    row[nc]=(-r0[i])%p if i in F16 else 0
    if any(row): A.append(row); ridx.append(i)
print(f"nonzero rows: {len(A)} ({time.time()-t0:.0f}s), eliminating...",flush=True)
piv=0
for c in range(nc):
    sel=-1
    for r in range(piv,len(A)):
        if A[r][c]%p: sel=r;break
    if sel<0: continue
    A[piv],A[sel]=A[sel],A[piv]; ridx[piv],ridx[sel]=ridx[sel],ridx[piv]
    ivv=inv(A[piv][c]); A[piv]=[(x*ivv)%p for x in A[piv]]
    for r in range(len(A)):
        if r!=piv and A[r][c]%p:
            f=A[r][c]; A[r]=[(A[r][k]-f*A[piv][k])%p for k in range(nc+1)]
    piv+=1
incons=[ridx[r] for r in range(len(A)) if all(A[r][k]%p==0 for k in range(nc)) and A[r][nc]%p]
print(f"rank={piv}, INCONSISTENT rows={len(incons)} ({time.time()-t0:.0f}s)")
if not incons: print("*** MOD-P CONSISTENT over this pool - witness reachable mod p, obstruction is z-lift ***")
else: print(f"mod-p INCONSISTENT even over {len(poolL)} frees. sample: {incons[:6]}")
