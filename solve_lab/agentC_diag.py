#!/usr/bin/env python3
"""Diagnose the breakage after applying a null-space (da,dc) move: are broken eqs mod-p or ℤ-only?
in the constraint set? linear or nonlinear in handles? core status."""
import json, sys
from agentC_common import (p, gates, order, definer, gcode, forward, val, freeinp, ns, lines,
                           eqcode, eqvars, load_best, CORE, NVARS, pinned, rootcode_of, inv,
                           is_qr, C1, C2, downstream_ks, partial_forward)
import agentC_poly as P
from collections import defaultdict

which = int(sys.argv[1]) if len(sys.argv)>1 else 0   # which cubic root
best=load_best(); forward()
V=json.load(open('agentC_Vdata.json'))
H=V['H']; Hidx={h:i for i,h in enumerate(H)}; nullbasis=V['nullbasis']
Hset=set(H)
x29322=val[29322]%p; x1326=val[1326]%p; x27713=val[27713]%p; x33469=val[33469]%p; x3558=val[3558]%p
T0=val[6671]%p
den=(x29322+x1326)%p
beta=((x27713+x3558)*inv(den))%p; alpha=((-T0)*inv(den))%p
f1=P.pmul([x33469,1], P.pmul([x29322,(-1)%p],[x29322,(-1)%p]))
lin=[(x3558-alpha)%p,(-beta)%p]; f2=P.pmul(lin,lin); Spoly=P.psub(f1,f2)
roots=P.roots_mod_p(Spoly)
da=roots[which]; dc=(alpha+beta*da)%p
print(f"root {which}: da={da}\n dc={dc}")

# constraint set from closure (recompute quickly: eqs touched by H, currently satisfied)
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
cons=set()
for h in H:
    aff=set(eqbyvar.get(h,()))
    for k in downstream_ks(h): aff|=eqbyvar.get(order[k],set())
    cons|=(aff-F0)
print(f"|cons|={len(cons)}")

# baseline roots mod p of all sat eqs
base_root={i:eval(rootcode_of(i),ns)%p for i in range(len(lines))}
# apply move
for h in H:
    d=(da*nullbasis[0][Hidx[h]] + dc*nullbasis[1][Hidx[h]])%p
    val[h]=val[h]+d
forward(); ns['v']=val
# quotient handles
q_ok={'L1':val[11150]%p==0,'L3':(537773*val[37758])%p==0,'L2':val[25739]%(6672769*p)==0}
print(f"quotient prereqs (L*≡0 mod p): {q_ok}")
if val[11150]%p==0: val[30317]=-(val[11150]//p)
if (537773*val[37758])%p==0: val[2936]=(537773*val[37758])//p
if val[25739]%(6672769*p)==0: val[5146]=val[25739]//(6672769*p)
forward(); ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
core=[i for i in F if i in CORE]; nc=[i for i in F if i not in CORE]
print(f"satisfied {len(lines)-len(F)}/{len(lines)}; core-fail={len(core)}; noncore-fail={len(nc)}: {sorted(nc)}")
print(f"S%p={val[35389]%p}, T%p={val[6671]%p}")

# diagnose each noncore break
print("\n-- noncore break diagnosis --")
for i in sorted(nc):
    rp=eval(rootcode_of(i),ns)%p           # root mod p now
    rz=eval(rootcode_of(i),ns)              # root over Z
    incons = i in cons
    modp_broke = rp!=0
    z_only = (rp==0 and rz!=0)
    print(f"  eq {i}: in_cons={incons}, modp_root={'NONZERO' if modp_broke else 0}, "
          f"z_only_carry={z_only}, |root|~{len(str(abs(rz)))}digits, "
          f"div_by_p={rz%p==0}")
# core diagnosis: how many core are 'downstream' of noncore breaks vs independent
print("\n-- core status: which core eqs fail --")
print(f"core fail: {sorted(core)}")
