#!/usr/bin/env python3
"""Grow unchecked-free pool, adding compensators for inconsistent eqs. Track consistency trend."""
import heal_harness as H
import pickle,time
from collections import defaultdict
p=H.p
CK=pickle.load(open('checked.pkl','rb')); checked=CK['checked']
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
ns={'v':H.val,'__builtins__':{}}
def inv(a): return pow(a%p,p-2,p)
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
base=H.val[:]
def jaccol(w,Et,r0):
    H.val[w]+=1
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
    out={}
    for i in Et:
        d=(eval(H.eqcode[i],ns)-r0[i])%p
        if d: out[i]=d
    H.val[w]=base[w]
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
    return out
# free-ancestors per eq (unchecked only)
eq_ufree=[]
for i in range(len(H.eqcode)):
    s=set()
    for w in H.eqvars[i]:
        if w in H.freeinp and w not in checked: s.add(w)
        s|=(H.anc.get(w,set()))
    eq_ufree.append((s&H.freeinp)-checked)
F16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
pool=set()
for i in F16: pool|=eq_ufree[i]
t0=time.time()
for git in range(8):
    poolL=sorted(pool)
    descset=set()
    for w in poolL: descset.update(H.order[k] for k in desc_of[w])
    descset|=set(poolL)
    Et=[i for i in range(len(H.eqcode)) if H.eqvars[i]&descset]
    r0={i:eval(H.eqcode[i],ns) for i in Et}
    col={w:j for j,w in enumerate(poolL)}; nc=len(poolL)
    Jr=defaultdict(dict)
    for w in poolL:
        for i,c in jaccol(w,Et,r0).items(): Jr[i][col[w]]=c
    rows=[(i,Jr[i],(-r0[i])%p) for i in Et if Jr[i] or r0[i]%p]
    A=[]; ridx=[]
    for i,d,rhs in rows:
        row=[0]*(nc+1)
        for j,c in d.items(): row[j]=c
        row[nc]=rhs; A.append(row); ridx.append(i)
    piv=0
    for c in range(nc):
        sel=-1
        for r in range(piv,len(A)):
            if A[r][c]%p: sel=r;break
        if sel<0: continue
        A[piv],A[sel]=A[sel],A[piv]; ridx[piv],ridx[sel]=ridx[sel],ridx[piv]
        iv=inv(A[piv][c]); A[piv]=[(x*iv)%p for x in A[piv]]
        for r in range(len(A)):
            if r!=piv and A[r][c]%p:
                f=A[r][c]; A[r]=[(A[r][k]-f*A[piv][k])%p for k in range(nc+1)]
        piv+=1
    inc=[ridx[r] for r in range(len(A)) if all(A[r][k]%p==0 for k in range(nc)) and A[r][nc]%p]
    print(f"grow {git}: pool={nc}, touched={len(Et)}, rank={piv}, inconsistent={len(inc)} ({time.time()-t0:.0f}s)")
    if not inc:
        print("CONSISTENT! (first-order) — unchecked heal feasible"); break
    add=set()
    for i in inc: add|=eq_ufree[i]
    if add<=pool:
        print("pool closed but inconsistent"); break
    pool|=add
    if len(pool)>800: print("pool cap reached"); break
