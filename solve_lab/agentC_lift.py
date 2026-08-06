#!/usr/bin/env python3
"""Solve S(da,dc)=T(da,dc)=0 mod p on the 2-dim achievable subspace, apply the null-space move,
set quotient handles, and measure breakage (mod p and exact ℤ)."""
import json, time
from agentC_common import (p, gates, order, definer, gcode, forward, val, freeinp, ns, lines,
                           eqcode, eqvars, load_best, CORE, NVARS, pinned, rootcode_of, inv,
                           is_qr, sqrt_mod, C1, C2)
import agentC_poly as P

best=load_best(); forward()
V=json.load(open('agentC_Vdata.json'))
H=V['H']; Hidx={h:i for i,h in enumerate(H)}
projbasis=V['projbasis']; nullbasis=V['nullbasis']
dimV=V['dimV']
assert dimV==2 and projbasis==[[1,0,0],[0,0,1]], (dimV,projbasis)
# current deep residues
x29322=val[29322]%p; x1326=val[1326]%p; x27713=val[27713]%p; x33469=val[33469]%p; x3558=val[3558]%p
S0=val[35389]%p; T0=val[6671]%p
print(f"current S={S0}, T={T0}")
# T(da,dc) = T0 - da*(x27713+x3558) + dc*(x29322+x1326) = 0  -> dc = alpha + beta*da
den=(x29322+x1326)%p
print(f"x29322+x1326 mod p = {den}")
assert den!=0, "degenerate: x29322+x1326==0, handle separately"
beta=((x27713+x3558)*inv(den))%p
alpha=((-T0)*inv(den))%p
# verify T linear form: pick random da, dc=alpha+beta*da => T==0
import random
for _ in range(3):
    da=random.randrange(p); dc=(alpha+beta*da)%p
    T=(( (x27713+dc)*((x29322-da)%p) - ((x3558-dc)%p)*((x1326+da)%p) ))%p
    assert T==0, T
print("T-elimination verified (dc = alpha + beta*da)")
# S(da) with dc=alpha+beta*da:  (x33469+da)*(x29322-da)^2 - (x3558-dc)^2
# build as polynomial in da (little-endian)
f1=P.pmul([x33469,1], P.pmul([x29322,(-1)%p],[x29322,(-1)%p]))   # (x33469+da)(x29322-da)^2
# x3558 - dc = (x3558-alpha) + (-beta)*da
lin=[(x3558-alpha)%p,(-beta)%p]
f2=P.pmul(lin,lin)
Spoly=P.psub(f1,f2)
print(f"S(da) cubic coeffs (deg {P.pdeg(Spoly)}): {Spoly}")
roots=P.roots_mod_p(Spoly)
print(f"cubic roots (da candidates): {len(roots)}")
for da in roots:
    dc=(alpha+beta*da)%p
    # verify S,T ==0
    A=(x33469+da)%p; X29=(x29322-da)%p; X35=(x3558-dc)%p; X13=(x1326+da)%p; X27=(x27713+dc)%p
    S=(A*X29*X29 - X35*X35)%p
    T=(X27*X29 - X35*X13)%p
    print(f"  da={da}\n  dc={dc}\n   -> S={S}, T={T}  (x_33469 QR now: {is_qr(A)})")

if not roots:
    print("NO cubic root -> try alternative param or regime; exiting")
    import sys; sys.exit(0)

# --- apply the FIRST root and measure breakage ---
def apply_and_measure(da,dc,label):
    # snapshot
    snap={h:val[h] for h in H}
    for h in H:
        d=(da*nullbasis[0][Hidx[h]] + dc*nullbasis[1][Hidx[h]])%p
        val[h]=(val[h]+d)
    forward(); ns['v']=val
    Snew=val[35389]%p; Tnew=val[6671]%p
    # set quotient handles
    if val[11150]%p==0: val[30317]=-(val[11150]//p)
    if (537773*val[37758])%p==0: val[2936]=(537773*val[37758])//p
    if val[25739]%(6672769*p)==0: val[5146]=val[25739]//(6672769*p)
    forward(); ns['v']=val
    F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    core=[i for i in F if i in CORE]; nc=[i for i in F if i not in CORE]
    print(f"[{label}] da,dc applied: S%p={Snew}, T%p={Tnew}")
    print(f"   satisfied {len(lines)-len(F)}/{len(lines)}; core-fail={len(core)}; noncore-fail={len(nc)}")
    print(f"   noncore (first 40): {sorted(nc)[:40]}")
    print(f"   L2/p % 6672769 = {(val[25739]//p)%6672769}")
    # restore
    for h in H: val[h]=snap[h]
    forward()
    return len(nc), core, nc

da0=roots[0]; dc0=(alpha+beta*da0)%p
apply_and_measure(da0,dc0,"root0")
