#!/usr/bin/env python3
"""Newton heal of the 16 ripple using ONLY unchecked frees (preserve equality checks).
Keep all touched equations satisfied. Iterate for nonlinearity."""
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
# descendants per free (incremental)
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
def recompute(ws):
    ks=set()
    for w in ws: ks.update(desc_of[w])
    for k in sorted(ks): H.val[H.order[k]]=eval(H.gcode[k],ns)
# pool: unchecked frees with leverage (compute once) - use a generous set
F16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
eq_anc=set()
for i in F16:
    for w in H.eqvars[i]:
        if w in H.freeinp: eq_anc.add(w)
        eq_anc|=H.anc.get(w,set())
pool=sorted((eq_anc&H.freeinp)-checked)
print(f"unchecked pool: {len(pool)}")
allfails=lambda: [i for i,c in enumerate(H.eqcode) if eval(c,ns)!=0]
F=allfails(); print(f"start fails: {len(F)}")
for it in range(15):
    # touched eqs = eqs whose vars intersect pool-descendants
    descset=set()
    for w in pool: descset.update(H.order[k] for k in desc_of[w])
    descset|=set(pool)
    Et=[i for i in range(len(H.eqcode)) if H.eqvars[i]&descset]
    r0={i:eval(H.eqcode[i],ns) for i in Et}
    # Jacobian
    col={w:j for j,w in enumerate(pool)}; nc=len(pool)
    Jrows=defaultdict(dict)
    base=H.val[:]
    for w in pool:
        H.val[w]+=1; recompute([w])
        for i in Et:
            d=(eval(H.eqcode[i],ns)-r0[i])%p
            if d: Jrows[i][col[w]]=d
        H.val[:]=base[:]; recompute([w])
    # build system: Et rows, target -r0 (for all Et, since we want all = 0)
    rows=[]
    for i in Et:
        if Jrows[i] or (r0[i]%p): rows.append((i,Jrows[i],(-r0[i])%p))
    A=[]; ridx=[]
    for i,d,rhs in rows:
        row=[0]*(nc+1)
        for j,c in d.items(): row[j]=c
        row[nc]=rhs; A.append(row); ridx.append(i)
    piv=0; where=[-1]*nc
    for c in range(nc):
        sel=-1
        for r in range(piv,len(A)):
            if A[r][c]%p: sel=r;break
        if sel<0: continue
        A[piv],A[sel]=A[sel],A[piv]
        iv=inv(A[piv][c]); A[piv]=[(x*iv)%p for x in A[piv]]
        for r in range(len(A)):
            if r!=piv and A[r][c]%p:
                f=A[r][c]; A[r]=[(A[r][k]-f*A[piv][k])%p for k in range(nc+1)]
        where[c]=piv; piv+=1
    incons=sum(1 for r in range(len(A)) if all(A[r][k]%p==0 for k in range(nc)) and A[r][nc]%p)
    if incons:
        print(f"it{it}: rank={piv}, INCONSISTENT rows={incons} -> unchecked heal blocked"); break
    d=[0]*nc
    for c in range(nc):
        if where[c]>=0: d[c]=A[where[c]][nc]%p
    # apply (lift to signed small)
    for w in pool:
        dv=d[col[w]]
        if dv>p//2: dv-=p
        H.val[w]=(H.val[w]+dv)
    recompute(pool)
    F=allfails()
    print(f"it{it}: rank={piv}, applied, now {len(F)} fails")
    if len(F)==0: print("SOLVED!!!"); import json; json.dump({f'x_{i}':H.val[i] for i in range(H.NVARS)},open('SOLVED_full.json','w')); break
