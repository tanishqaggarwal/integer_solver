#!/usr/bin/env python3
"""Extract the locked coupling ratio r1 = (dG2)/(dG1) on null(J_sat), and verify the
mod-p obstruction certificate:  G2 != r1*G1 (mod p)  => cannot close both gaps."""
import pickle,os,sys
from ge import solve_prime,inv
sp=os.path.dirname(os.path.abspath(__file__))
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
d=pickle.load(open(os.path.join(sp,f'jac_{p}.pkl'),'rb'))
rows=d['rows']; g1=d['g1']; g2=d['g2']
G1=d['G1i']%p; G2=d['G2i']%p
g1grad=g1[0]; g2grad=g2[0]
ZERO=({},0)
def dot(row,sol): return sum(c*sol.get(k,0) for k,c in row.items())%p

# Solve {J_sat delta=0, gradG1.delta = 1}, then r1 = gradG2.delta
cons,sol=solve_prime(rows, (g1grad,1), ZERO, p, d['ncol'], want_sol=True, verbose=False)
assert cons, "gap1(rhs=1) should be consistent"
chk1=dot(g1grad,sol)
r1=dot(g2grad,sol)
print(f"gradG1.delta = {chk1} (want 1)")
print(f"r1 = gradG2.delta on null(J_sat) with gradG1.delta=1:")
print(f"  r1 = {r1}")

# cross-check by driving gradG1.delta = 5 -> gradG2.delta should be 5*r1
cons2,sol2=solve_prime(rows, (g1grad,5), ZERO, p, d['ncol'], want_sol=True, verbose=False)
print(f"  cross-check drive=5: gradG2.delta={dot(g2grad,sol2)}  5*r1={(5*r1)%p}  match={dot(g2grad,sol2)==(5*r1)%p}")

# The obstruction: to close, need gradG1.delta=-G1 => gradG2.delta = -r1*G1, but need -G2.
lhs=(r1*G1)%p          # forced G2-change magnitude if we set gradG1.delta=-G1 (sign folded)
print(f"\nG1 mod p = {G1}")
print(f"G2 mod p = {G2}")
print(f"r1*G1 mod p = {lhs}")
obstruction=(G2 - lhs)%p
print(f"OBSTRUCTION  Q = (G2 - r1*G1) mod p = {obstruction}")
print(f"=> closing both gaps requires Q==0.  Q {'==0 (would be FEASIBLE)' if obstruction==0 else '!=0  => INFEASIBLE mod p'}")

# Also compute r2 (drive gradG2, read gradG1) as sanity: should be inverse relation r1*r2==1 if both nonzero
cons3,sol3=solve_prime(rows, (g2grad,1), ZERO, p, d['ncol'], want_sol=True, verbose=False)
if cons3:
    r2=dot(g1grad,sol3)
    print(f"\nr2 = gradG1.delta on null(J_sat) with gradG2.delta=1: {r2}")
    print(f"  r1*r2 mod p = {(r1*r2)%p}  (==1 confirms exact 1-dim locked coupling)")
