#!/usr/bin/env python3
import pickle
p=2**256-2**32-977
D=pickle.load(open('heal16.pkl','rb'))
pool=D['pool']; Etouch=D['Etouch']; r0=D['r0']; J=D['J']; F16=set(D['F16'])
col={w:j for j,w in enumerate(pool)}; nc=len(pool)
# build rows: for each eq in Etouch, coeffs over pool, rhs = target
# target: F16 -> want new resid 0 => J*d = -r0 ; others -> want stay 0 => J*d = 0
rows=[]
for i in Etouch:
    row=[0]*nc
    any_=False
    for w in pool:
        c=J.get((i,w),0)
        if c: row[col[w]]=c%p; any_=True
    rhs=(-r0[i])%p if i in F16 else 0
    if any_ or rhs: rows.append((row,rhs,i))
print(f"rows with content: {len(rows)}, unknowns: {nc}")
# Gaussian elimination mod p on augmented [row|rhs]
A=[r[0][:]+[r[1]] for r in rows]
idxs=[r[2] for r in rows]
def inv(a): return pow(a%p,p-2,p)
piv_r=0; where=[-1]*nc
for cljj in range(nc):
    sel=-1
    for r in range(piv_r,len(A)):
        if A[r][cljj]%p!=0: sel=r; break
    if sel<0: continue
    A[piv_r],A[sel]=A[sel],A[piv_r]
    iv=inv(A[piv_r][cljj])
    A[piv_r]=[(x*iv)%p for x in A[piv_r]]
    for r in range(len(A)):
        if r!=piv_r and A[r][cljj]%p!=0:
            f=A[r][cljj]
            A[r]=[(A[r][k]-f*A[piv_r][k])%p for k in range(nc+1)]
    where[cljj]=piv_r; piv_r+=1
# consistency: any row all-zero coeff but nonzero rhs?
incons=[]
for r in range(len(A)):
    if all(A[r][k]%p==0 for k in range(nc)) and A[r][nc]%p!=0:
        incons.append(idxs[r] if r<len(idxs) else '?')
print(f"rank={piv_r}, inconsistent rows: {len(incons)}")
if incons:
    print("INCONSISTENT — some kept-eqs conflict. sample bad eq idx:", incons[:10])
else:
    print("CONSISTENT! solving...")
    d=[0]*nc
    for cljj in range(nc):
        if where[cljj]>=0: d[cljj]=A[where[cljj]][nc]%p
    sol={pool[j]:d[j] for j in range(nc) if d[j]!=0}
    print(f"solution nonzero knobs: {len(sol)}")
    pickle.dump(sol, open('heal16_sol.pkl','wb'))
    print("saved heal16_sol.pkl:", {k:'..' for k in list(sol)[:8]})
