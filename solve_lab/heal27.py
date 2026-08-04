#!/usr/bin/env python3
"""Heal the 27 core-ripple equations using FINE-GRAINED free slacks (free vars that multiply
1-valued or small non-p vars). These can inject arbitrary sub-p corrections."""
import heal_harness as H
import pickle
from collections import defaultdict
p=H.p
def inv(a): return pow(a%p,p-2,p)
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=v013.get(v,0)
H.val[14853]=vA[14853]; H.val[31339]=vA[31339]
H.forward()
ns={'v':H.val,'__builtins__':{}}
F27=[i for i in range(len(H.eqcode)) if eval(H.eqcode[i],ns)!=0]
print(f"core-only: {len(F27)} fails")
# descendants per free
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
# Pool: free-zero vars in the 27's cone that are NOT wire and NOT the p-granular gap slacks... 
# actually include everything free that affects the 27, EXCEPT x_4432,x_7068 (protect G1,G2) and core knobs
PROTECT={4432,7068,14853,31339,12186,16742}
pool=set()
for i in F27:
    for w in H.eqvars[i]:
        if w in H.freeinp: pool.add(w)
        pool|=(H.anc.get(w,set())&H.freeinp)
pool-=PROTECT
pool=sorted(pool)
print(f"pool (free, excl protect): {len(pool)}")
# Jacobian of the 27 wrt pool
base=H.val[:]
r0={i:eval(H.eqcode[i],ns)%p for i in F27}
col={w:j for j,w in enumerate(pool)}; nc=len(pool)
J=defaultdict(dict)
for w in pool:
    H.val[w]+=1
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
    for i in F27:
        d=(eval(H.eqcode[i],ns)-r0[i])%p
        if d: J[i][col[w]]=d
    H.val[w]=base[w]
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
# also must keep currently-satisfied eqs touched by pool satisfied
descset=set()
for w in pool: descset.update(H.order[k] for k in desc_of[w])
descset|=set(pool)
Et=[i for i in range(len(H.eqcode)) if i not in set(F27) and (H.eqvars[i]&descset)]
r0e={i:eval(H.eqcode[i],ns)%p for i in Et}
for w in pool:
    H.val[w]+=1
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
    for i in Et:
        d=(eval(H.eqcode[i],ns)-r0e[i])%p
        if d: J[i][col[w]]=d
    H.val[w]=base[w]
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],ns)
# system: F27 -> -r0 ; Et -> 0
rows=[]
allrhs={**{i:(-r0[i])%p for i in F27},**{i:0 for i in Et}}
for i in list(F27)+Et:
    if J[i] or (i in F27 and r0[i]):
        rows.append((i,J[i],allrhs[i]))
A=[];ridx=[]
for i,dd,rhs in rows:
    row=[0]*(nc+1)
    for j,c in dd.items(): row[j]=c
    row[nc]=rhs%p; A.append(row); ridx.append(i)
piv=0;where=[-1]*nc
for c in range(nc):
    sel=-1
    for r in range(piv,len(A)):
        if A[r][c]%p: sel=r;break
    if sel<0: continue
    A[piv],A[sel]=A[sel],A[piv];ridx[piv],ridx[sel]=ridx[sel],ridx[piv]
    ivv=inv(A[piv][c]);A[piv]=[(x*ivv)%p for x in A[piv]]
    for r in range(len(A)):
        if r!=piv and A[r][c]%p:
            f=A[r][c];A[r]=[(A[r][k]-f*A[piv][k])%p for k in range(nc+1)]
    where[c]=piv;piv+=1
incons=[ridx[r] for r in range(len(A)) if all(A[r][k]%p==0 for k in range(nc)) and A[r][nc]%p]
print(f"pool={nc}, touched-satisfied={len(Et)}, rank={piv}, INCONSISTENT={len(incons)}")
inc_in_27=[i for i in incons if i in set(F27)]
print(f"  inconsistent that are among the 27: {len(inc_in_27)} -> {inc_in_27}")
print(f"  inconsistent (satisfied eqs that break): {len([i for i in incons if i not in set(F27)])}")
