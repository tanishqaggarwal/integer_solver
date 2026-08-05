"""Exact GF(p) Jacobian of the 27 residual eqs (11 fails + 16 ripple) w.r.t. ALL free inputs.
Test first-order consistency of J*delta = -r mod p (can we close all 27 to first order?)."""
import heal_harness as H
from jac_lib import D, freelist, freeidx, NF
import json
p=H.p
FAILS=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIP=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
EQS=FAILS+RIP
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
# seed duals on ALL free inputs
vd=[None]*H.NVARS
for j in H.freeinp:
    vd[j]=D(H.val[j],{freeidx[j]:1})
ns={'v':vd,'__builtins__':{}}
for k,t in enumerate(H.order):
    r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
# build rows for the 27 eqs
rows=[]; rhs=[]
for i in EQS:
    rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
    if isinstance(rr,D): rv=rr.v; g=rr.g
    else: rv=rr%p; g={}
    rows.append(dict(g)); rhs.append((-rv)%p)
nr=len(rows)
# columns actually used
cols=sorted(set(c for g in rows for c in g))
print(f'{nr} eqs, {len(cols)} active free-input columns')
# residuals: how many nonzero (mod p)?
nzr=[i for i,b in enumerate(rhs) if b!=0]
print(f'nonzero residuals (mod p) among 27: {len(nzr)} -> eqs {[EQS[i] for i in nzr]}')
# Gaussian elimination mod p on augmented [J | -r], detect inconsistency
def inv(a): return pow(a%p,p-2,p)
R=[dict(g) for g in rows]; B=list(rhs)
colset=cols
rank=0; pivcols=[]
usedrow=[False]*nr
for c in colset:
    # find pivot row
    prow=-1
    for r in range(nr):
        if not usedrow[r] and R[r].get(c,0)!=0: prow=r; break
    if prow<0: continue
    usedrow[prow]=True; pivcols.append(c); rank+=1
    iv=inv(R[prow][c])
    # normalize
    R[prow]={k:(v*iv)%p for k,v in R[prow].items()}; B[prow]=(B[prow]*iv)%p
    for r in range(nr):
        if r!=prow and R[r].get(c,0)!=0:
            f=R[r][c]
            for k,v in R[prow].items():
                nv=(R[r].get(k,0)-f*v)%p
                if nv: R[r][k]=nv
                elif k in R[r]: del R[r][k]
            B[r]=(B[r]-f*B[prow])%p
# after elimination: inconsistent if some row has empty lhs but nonzero rhs
incon=[r for r in range(nr) if not R[r] and B[r]!=0]
print(f'rank(J)={rank}')
print(f'INCONSISTENT rows (0 = b!=0): {len(incon)}  -> {"FIRST-ORDER INFEASIBLE" if incon else "FIRST-ORDER FEASIBLE (mod p)"}')
if incon:
    # the inconsistency combination
    print('  wall confirmed: cannot close all 27 to first order mod p')
