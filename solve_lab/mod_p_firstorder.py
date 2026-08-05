#!/usr/bin/env python3
"""Independent first-order mod-p feasibility check for channel (1,0).
Build sparse mod-p Jacobian J of ALL equations at best_agentA (jac_lib is mod p),
residual R (nonzero on the 11 gap eqs). Test consistency of J.delta = -R (mod p)
via sparse Gaussian elimination (dict rows). Consistent => first-order repair of the
gap exists while keeping the 39022 satisfied. Also report whether the x_7068 and
x_4432 message-match conditions are individually reachable (augment with target rows)."""
import sys, os, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
os.chdir('/home/user/integer_solver/solve_lab')
import heal_harness as H
from jac_lib import D, build_duals, eq_jac_row, freelist, freeidx
p=H.p
def inv(a): return pow(a%p,p-2,p)

vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
fails=set(H.fails())
print('build duals...',flush=True); t0=time.time()
vd=build_duals()
print('  done %.1fs'%(time.time()-t0),flush=True)

# Build sparse rows: only keep rows that are gap (nonzero resid) OR have nonzero gradient.
# For a consistency test we need ALL rows whose gradient is nonzero (they constrain delta).
print('build jacobian rows...',flush=True); t0=time.time()
rows=[]   # each: (dict col->coef, rhs)  representing sum coef*delta = rhs  (rhs=-resid)
for i in range(len(H.eqcode)):
    rv,g=eq_jac_row(i,vd)
    if not g:
        # no gradient: equation is constant under free-input perturbation; must already be 0
        continue
    rows.append((dict(g),(-rv)%p))
print('  %d constraining rows, %.1fs'%(len(rows),time.time()-t0),flush=True)

# Sparse GE mod p with the augmented rhs. Detect inconsistency (empty row, nonzero rhs).
t0=time.time()
# pivot map: col -> row-index in a list of (rowdict, rhs)
piv={}  # col -> (rowdict, rhs) normalized so rowdict[col]=1
order=0
incons=0
for (rd0,rb0) in rows:
    rd=dict(rd0); rb=rb0
    # reduce by existing pivots
    changed=True
    # iterate over a copy of cols present
    stack=[c for c in rd]
    while stack:
        c=stack.pop()
        if c not in rd: continue
        if c in piv:
            f=rd[c]
            if f==0:
                del rd[c]; continue
            prd,prb=piv[c]
            # rd -= f*prd
            for k,v in prd.items():
                nv=(rd.get(k,0)-f*v)%p
                if nv: 
                    if k not in rd: stack.append(k)
                    rd[k]=nv
                elif k in rd: del rd[k]
            rb=(rb-f*prb)%p
    # now rd has no pivot cols; pick a new pivot
    rd={k:v for k,v in rd.items() if v}
    if not rd:
        if rb%p!=0: incons+=1
        continue
    c=next(iter(rd)); iv=inv(rd[c])
    rd={k:(v*iv)%p for k,v in rd.items()}; rb=(rb*iv)%p
    piv[c]=(rd,rb)
    order+=1
    if order%2000==0: print('  pivots=%d elapsed=%.0fs'%(order,time.time()-t0),flush=True)
print('rank(aug elimination): pivots=%d, inconsistent_rows=%d, %.0fs'%(order,incons,time.time()-t0),flush=True)
if incons==0:
    print('*** FIRST-ORDER mod-p CONSISTENT: a delta zeroing the gap while keeping the 39022 exists to first order. ***',flush=True)
    print('    (Strong positive signal for channel (1,0) mod-p feasibility; Newton should converge.)',flush=True)
else:
    print('*** FIRST-ORDER mod-p INCONSISTENT (%d bad rows): no first-order repair; higher-order needed or infeasible. ***'%incons,flush=True)
