#!/usr/bin/env python3
"""Linearized mod-p solve of the 39-equation bit system. Unknowns t_b=b, s_b=b^2, u_ab=a*b.
Solve linear system mod p; check consistency s_b==t_b^2, u_ab==t_a*t_b."""
import pickle
p=2**256-2**32-977
D=pickle.load(open('bitpolys.pkl','rb'))
polys=D['polys']; BE=D['BE']; BITS=D['BITS']
# unknown index
tvar={b:('t',b) for b in BITS}
svar={b:('s',b) for b in BITS}
uvars={}
for i in BE:
    for (a,b) in polys[i]['cross']:
        uvars[(min(a,b),max(a,b))]=('u',a,b)
cols=[]
for b in BITS: cols.append(('t',b))
for b in BITS: cols.append(('s',b))
for k in uvars: cols.append(('u',)+k)
cidx={c:j for j,c in enumerate(cols)}; nc=len(cols)
print(f"unknowns: {nc} (13 t + 13 s + {len(uvars)} u)")
# build rows
rows=[]
for i in BE:
    P=polys[i]; row=[0]*nc; 
    for b,l in P['lin'].items():
        if l%p: row[cidx[('t',b)]]=(row[cidx[('t',b)]]+l)%p
    for b,q in P['quad'].items():
        if q%p: row[cidx[('s',b)]]=(row[cidx[('s',b)]]+q)%p
    for (a,b),cr in P['cross'].items():
        row[cidx[('u',min(a,b),max(a,b))]]=(cr)%p
    rhs=(-P['c0'])%p
    rows.append((row,rhs,i))
# gaussian elim mod p
def inv(a): return pow(a%p,p-2,p)
A=[r[0][:]+[r[1]] for r in rows]; ridx=[r[2] for r in rows]
piv=0; where=[-1]*nc
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
    where[c]=piv; piv+=1
incons=[ridx[r] for r in range(len(A)) if all(A[r][k]%p==0 for k in range(nc)) and A[r][nc]%p]
print(f"rank={piv}, free unknowns={nc-piv}, inconsistent rows={len(incons)}")
if incons:
    print("mod-p linear system INCONSISTENT:",incons[:8])
else:
    # extract solution (set free vars=0)
    sol=[0]*nc
    for c in range(nc):
        if where[c]>=0: sol[c]=A[where[c]][nc]%p
    tb={b:sol[cidx[('t',b)]] for b in BITS}
    sb={b:sol[cidx[('s',b)]] for b in BITS}
    print("free unknowns present -> particular solution (free=0). Checking s_b==t_b^2 mod p:")
    ok=0
    for b in BITS:
        good=(sb[b]==(tb[b]*tb[b])%p)
        ok+=good
        print(f"  bit x_{b}: t={'0' if tb[b]==0 else 'set'}, s==t^2 mod p: {good}")
    print(f"consistency s==t^2: {ok}/{len(BITS)}")
    pickle.dump({'tb':tb,'sb':sb,'where':where,'A':A,'cols':cols,'cidx':cidx,'nc':nc}, open('bitsol.pkl','wb'))
