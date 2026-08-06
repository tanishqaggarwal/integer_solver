#!/usr/bin/env python3
"""Beyond-first-order confirmation over ALL 8583 free inputs.
 Phase A: mod-p Newton driving {39022 sat eqs =0, G1=0} to an EXACT mod-p solution from best_agentA.
 Then report G2 at that solution and test whether G2 can be closed there (keeping sat & G1)."""
import os,sys
os.environ.setdefault("JAC_M","115792089237316195423570985008687907853269984665640564039457584007908834671663")
import heal_harness as H
import jac_lib_m as J
from ge import solve_prime
p=J.p
GAP=set([2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125])
vA=H.loadd('best_agentA_39022.json')
fv={j:vA.get(j,0)%p for j in H.freeinp}
ZERO=({},0)
def relin(fv):
    for j in H.freeinp: H.val[j]=fv[j]
    vd=J.build_duals()
    gd1=vd[7068]-vd[2099]-7376877*vd[642]
    gd2=vd[4432]-vd[19964]-vd[28730]
    rows=[]; nsatbad=0
    for i in range(len(H.eqcode)):
        if i in GAP: continue
        r=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
        if isinstance(r,J.D):
            if r.v%p: nsatbad+=1
            g={k:c%p for k,c in r.g.items() if c%p}
            if g or r.v%p: rows.append((g,(-r.v)%p))
    return vd,gd1,gd2,rows,nsatbad
def apply(fv,sol):
    for j in H.freeinp:
        c=J.freeidx[j]
        if c in sol: fv[j]=(fv[j]+sol[c])%p
print("Phase A: Newton on {39022 sat eqs, G1=0} (a CONSISTENT system):",flush=True)
G2v=None
for it in range(8):
    vd,gd1,gd2,rows,nbad=relin(fv)
    G1v=gd1.v%p; G2v=gd2.v%p
    print(f" iter {it}: sat-violations={nbad:6d}  G1!=0={G1v!=0}  G2!=0={G2v!=0}  G2[:14]={str(G2v)[:14]}",flush=True)
    if nbad==0 and G1v==0:
        print("  -> reached EXACT mod-p solution of {39022 sat, G1=0}",flush=True); break
    g1row=({k:c%p for k,c in gd1.g.items() if c%p},(-G1v)%p)
    cons,sol=solve_prime(rows,g1row,ZERO,p,J.NF,want_sol=True,verbose=False)
    if not cons:
        print("  Newton step INCONSISTENT -- unexpected"); break
    apply(fv,sol)
print(f"\nAt the {{39022 sat, G1=0}} mod-p solution:  G2 mod p = {G2v}",flush=True)
print(f"  G2==0 ? {G2v==0}  ->  {'BOTH gaps closed (FEASIBLE mod p)' if G2v==0 else 'G2 pinned NONZERO: obstruction persists beyond first order'}",flush=True)
# At this new point, can we close G2 keeping {sat, G1=0} ?
vd,gd1,gd2,rows,nbad=relin(fv)
g1keep=({k:c%p for k,c in gd1.g.items() if c%p},0)
g2close=({k:c%p for k,c in gd2.g.items() if c%p},(-gd2.v)%p)
rows2=list(rows)+[g1keep]
cons,_=solve_prime(rows2,g2close,ZERO,p,J.NF,want_sol=False,verbose=False)
print(f"  Can close G2 at this new solution (keeping sat & G1=0)? {'YES' if cons else 'NO -> INCONSISTENT again (obstruction is not a first-order artifact of best_agentA)'}",flush=True)
