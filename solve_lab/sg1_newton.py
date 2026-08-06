"""Damped / Levenberg-Marquardt Newton on core residual (x_29322 mod p, x_3558 mod p) over
unchecked frees, with line-search on total integer fail count. Tests whether damping converges
where Agent B's undamped Newton diverged."""
import sys, json, time, math
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import heal_harness as H
p = H.p; val = H.val
from collections import defaultdict
ns = {'v': val, '__builtins__': {}}
eqcode=H.eqcode; eqvars=H.eqvars

efs=[]
for i in range(len(eqvars)):
    s=set()
    for var in eqvars[i]:
        if var in H.freeinp: s.add(var)
        else: s|=H.anc.get(var,set())
    efs.append(s)
sm=json.load(open('sg1_slavemap.json')); checked=set(int(v) for v in sm)

def load(path):
    d=H.loadd(path)
    for v in H.freeinp: val[v]=d.get(v,0)
    H.forward()
def nfails():
    return sum(1 for c in eqcode if eval(c,ns)!=0)
def resid():
    return (val[29322]%p, val[3558]%p)

load('best/new_instance_partial_39013.json')
r0=resid(); nf0=nfails()
print(f"39013: resid=(x29%p,x35%p)=({r0[0]!=0},{r0[1]!=0}) nfails={nf0}")

# knob set: unchecked frees that (probed) move x_29322 or x_3558 mod p
# use the slave-cone candidates found earlier + broaden
cand = [112,2055,2527,3114,5669,7151,14485,20454,24134,26277,27738,30163,33612,33787,35537,5910,24365]
cand = [f for f in cand if f in H.freeinp and f not in checked]

def jac_modp(knobs):
    """finite-diff Jacobian of (x29%p,x35%p) wrt each knob (bump +1)."""
    base=resid()
    J=[]
    for f in knobs:
        old=val[f]; val[f]=old+1; H.forward()
        r=resid(); val[f]=old; H.forward()
        J.append(((r[0]-base[0])%p,(r[1]-base[1])%p))
    return J,base

J,base=jac_modp(cand)
print(f"\nmod-p Jacobian rows (knob: d(x29%p), d(x35%p)):")
nz29=sum(1 for a,b in J if a!=0); nz35=sum(1 for a,b in J if b!=0)
print(f"  {nz29}/{len(cand)} knobs move x29%p; {nz35}/{len(cand)} move x35%p")

# Try a damped Newton step: pick knob with largest |d(x29%p)| effect, compute the EXACT integer
# bump k to zero x29%p: k ≡ -base0 * inv(d) mod p. Report its magnitude (the CVP step size).
import math
for f,(d29,d35) in zip(cand,J):
    if d29!=0:
        k=(-base[0]*pow(d29,p-2,p))%p
        # smallest signed representative
        ks = k if k< p//2 else k-p
        print(f"  knob {f}: exact bump to zero x29%p = {k}  (|signed|~{abs(ks).bit_length()} bits)")
        break

# Damped line-search: apply fractional steps of that bump, measure fail count growth
print(f"\nDamped homotopy: push knob {f} by fractions of the zeroing bump, watch fails:")
old=val[f]
for frac in [1e-40,1e-30,1e-20,1e-10,1e-3,0.01,0.1,0.5,1.0]:
    step=int(k*frac)
    val[f]=old+step; H.forward()
    print(f"  frac={frac:.0e} step_bits={step.bit_length()}: x29%p_zero={val[29322]%p==0} nfails={nfails()}")
val[f]=old; H.forward()
