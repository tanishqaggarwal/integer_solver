#!/usr/bin/env python3
"""Explicit mod-p rank / consistency numbers."""
import pickle,os
from ge import solve_prime
sp=os.path.dirname(os.path.abspath(__file__))
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
d=pickle.load(open(os.path.join(sp,f'jac_{p}.pkl'),'rb'))
rows=d['rows']; g1=d['g1']; g2=d['g2']; ncol=d['ncol']
G1=d['G1i']%p; G2=d['G2i']%p
ZERO=({},0)
g1h=(g1[0],0); g2h=(g2[0],0)   # homogeneous versions for pure coefficient rank
def rank_of(rws,extra):
    # homogeneous: count pivots
    cons=[True]
    piv=[0]
    # reuse solve_prime pivot count by giving rhs 0 everywhere; rank = pivots
    # solve_prime returns consistent + sol; we need pivot count -> re-implement quick
    from ge import inv
    pivots={}; order=[]
    allr=list(rws)+list(extra)
    for g,_ in allr:
        r=dict(g);
        while True:
            P=[c for c in r if c in pivots]
            if not P: break
            c=min(P); f=r[c]; pr=pivots[c]
            for cc,vv in pr.items():
                nv=(r.get(cc,0)-f*vv)%p
                if nv: r[cc]=nv
                elif cc in r: del r[cc]
        if r:
            c=min(r); ic=inv(r[c],p)
            pivots[c]={cc:(vv*ic)%p for cc,vv in r.items()}; order.append(c)
    return len(order)
r_sat   = rank_of(rows, [])
r_sat_1 = rank_of(rows, [g1h])
r_sat_12= rank_of(rows, [g1h,g2h])
print(f"rank(J_sat)                 = {r_sat}")
print(f"rank([J_sat; gradG1])       = {r_sat_1}   (+{r_sat_1-r_sat}: gradG1 {'INDEPENDENT of' if r_sat_1>r_sat else 'IN'} rowspace(J_sat))")
print(f"rank([J_sat; gradG1; gradG2])= {r_sat_12}   (+{r_sat_12-r_sat_1}: gradG2 {'INDEPENDENT of' if r_sat_12>r_sat_1 else 'IN'} rowspace([J_sat;gradG1]))")
# inhomogeneous consistency
cons,_=solve_prime(rows,g1,g2,p,ncol,want_sol=False,verbose=False)
print(f"Inhomogeneous system [J_sat.δ=0, gradG1.δ=-G1, gradG2.δ=-G2] consistent? {cons}")
print(f"=> rank(coeff)={r_sat_12}, rank([coeff|rhs])={r_sat_12+(0 if cons else 1)}  (Rouche-Capelli: {'consistent' if cons else 'INCONSISTENT'})")
print(f"nullity(J_sat) = ncol - rank = {ncol} - {r_sat} = {ncol-r_sat}")
