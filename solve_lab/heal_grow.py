#!/usr/bin/env python3
"""Incremental forward-reconstruct + iterative pool growth heal of the 16 ripple eqs.
Pin x_4432,x_7068 (protect 11) and core knobs. Grow compensator pool until consistent or capped."""
import heal_harness as H
import time
from collections import defaultdict
p=H.p
# precompute descendants per free input (topo-ordered target list)
t0=time.time()
topo_index={t:k for k,t in enumerate(H.order)}
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]:
        desc_of[w].append(k)   # k is topo index into H.order/gcode
print(f"precomputed desc_of in {time.time()-t0:.1f}s")
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F16=set(H.fails())
print(f"start fails: {len(F16)}")
PIN={4432,7068,14853,31339}
ns={'v':H.val,'__builtins__':{}}
gcode=[H.gcode[k] for k in range(len(H.order))]
order=H.order
def bump_eval(w, targets_needed):
    """bump w by 1, re-eval descendants, return {eq: new_resid} for eqs in targets_needed, then restore."""
    H.val[w]+=1
    for k in desc_of[w]: H.val[order[k]]=eval(gcode[k],ns)
    out={i:eval(H.eqcode[i],ns) for i in targets_needed}
    H.val[w]-=1
    for k in desc_of[w]: H.val[order[k]]=eval(gcode[k],ns)
    return out
# eq -> free ancestors (for pool growth + touched detection)
eq_free=[]
for i in range(len(H.eqcode)):
    s=set()
    for wv in H.eqvars[i]:
        if wv in H.freeinp: s.add(wv)
        s|=H.anc.get(wv,set())
    eq_free.append(s)
free_eqs=defaultdict(set)
for i in range(len(H.eqcode)):
    for wv in eq_free[i]: free_eqs[wv].add(i)

def inv(a): return pow(a%p,p-2,p)
def solve_pool(pool):
    pool=sorted(pool)
    # touched eqs = union of free_eqs[w] for w in pool
    Et=set()
    for w in pool: Et|=free_eqs[w]
    Et=sorted(Et)
    r0={i:eval(H.eqcode[i],ns)%p for i in Et}
    col={w:j for j,w in enumerate(pool)}; nc=len(pool)
    # Jacobian
    Jrows=defaultdict(dict)  # i -> {j:coef}
    for w in pool:
        newr=bump_eval(w,Et)
        for i in Et:
            d=(newr[i]-r0[i]*1)  # r0 is mod p but newr full; recompute delta mod p
        # redo properly:
    # (recompute with mod)
    for w in pool:
        newr=bump_eval(w,Et)
        for i in Et:
            d=(newr[i]-eval_full(i))%p
            if d: Jrows[i][col[w]]=d
    # build augmented rows
    rows=[]
    for i in Et:
        if not Jrows[i] and (i in F16):
            rows.append((i,{}, (-r0[i])%p))  # unfixable target
        elif Jrows[i]:
            rhs=(-r0[i])%p if i in F16 else 0
            rows.append((i,Jrows[i],rhs))
    # gaussian
    A=[]; ridx=[]
    for i,d,rhs in rows:
        row=[0]*(nc+1)
        for j,c in d.items(): row[j]=c%p
        row[nc]=rhs%p
        A.append(row); ridx.append(i)
    piv=0; where=[-1]*nc
    for c in range(nc):
        sel=-1
        for r in range(piv,len(A)):
            if A[r][c]%p!=0: sel=r;break
        if sel<0: continue
        A[piv],A[sel]=A[sel],A[piv]; ridx[piv],ridx[sel]=ridx[sel],ridx[piv]
        iv=inv(A[piv][c]); A[piv]=[(x*iv)%p for x in A[piv]]
        for r in range(len(A)):
            if r!=piv and A[r][c]%p!=0:
                f=A[r][c]; A[r]=[(A[r][k]-f*A[piv][k])%p for k in range(nc+1)]
        where[c]=piv; piv+=1
    incons=[ridx[r] for r in range(len(A)) if all(A[r][k]%p==0 for k in range(nc)) and A[r][nc]%p!=0]
    return sorted(pool),Et,piv,incons,A,where,nc
def eval_full(i): return eval(H.eqcode[i],ns)
pool=set(w for i in F16 for w in eq_free[i])-PIN
for git in range(12):
    poolL,Et,rank,incons,A,where,nc=solve_pool(pool)
    print(f"grow {git}: pool={len(poolL)}, touched={len(Et)}, rank={rank}, inconsistent={len(incons)} ({time.time()-t0:.0f}s)")
    if not incons:
        print("CONSISTENT!"); break
    add=set()
    for i in incons[:80]:
        add|=(eq_free[i]-PIN)
    if add<=pool:
        print("pool stopped growing but still inconsistent -> structural wall"); break
    pool|=add
