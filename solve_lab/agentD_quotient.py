#!/usr/bin/env python3
"""Task item 4: for a config, drive S=x_35389, T=x_6671 to 0 mod p via residue fixes on the
control-difference free inputs (x_14853, x_16742, x_12186), apply quotient handles
(x_30317=-L1//p, x_2936=537773*L3//p, x_5146=L2//(6672769*p)), re-forward, and measure count.
Quantifies the multi-role ripple damage."""
import sys, json
import agentD_harness as H
C1,C2=H.C1,H.C2
p=H.p
quad=sys.argv[1] if len(sys.argv)>1 else '11'
a7=24601; a34=2081
if quad=='11': ov={a7:1,a34:1,30213:C2,22162:C1,24468:C1,18956:C2}
elif quad=='10': ov={a7:1,30213:0,22162:0,24468:C1,18956:C2}
elif quad=='01': ov={a34:1,30213:0,22162:0,24468:C1,18956:C2}
elif quad=='00': ov={30213:0,22162:0,24468:C1,18956:C2}
r=H.run_config(ov, want_val=True)
val=r['val']; H.val=val; H.ns['v']=val
def recheck():
    ns={'v':val,'__builtins__':{}}
    F=[i for i in range(len(H.lines)) if eval(H.eqcode[i],ns)!=0]
    return F
print(f"[{quad}] base: sat={r['satisfied']} S%p={val[35389]%p} T%p={val[6671]%p} x15298={val[15298]}")
print(f"  x_29322={val[29322]%p} (=x_14853-x_12186), x_3558={val[3558]%p} (=x_24908-x_16742)")
print(f"  L1=x_11150%p={val[11150]%p}, L2=x_25739%p={val[25739]%p}, L3=x_37758%p={val[37758]%p}")

# Residue fix: make x_29322 ≡0 and x_3558 ≡0 mod p by adjusting the free-input sides.
# x_29322 = x_14853 - x_12186 ; reduce x_14853 by (x_29322 % p)
# x_3558  = x_24908 - x_16742 ; increase x_16742 by (x_3558 % p)  [x_16742 free]
d29322 = val[29322] % p
d3558  = val[3558] % p
print(f"\nApplying residue fix: x_14853 -= {d29322!=0}, x_16742 += {d3558!=0}")
val[14853] = (val[14853] - d29322)
val[16742] = (val[16742] + d3558)
H.forward();
print(f"  after fix (pre-quotient): S%p={val[35389]%p} T%p={val[6671]%p}")
print(f"  x_29322%p={val[29322]%p} x_3558%p={val[3558]%p}")
L1=val[11150]; L2=val[25739]; L3=val[37758]
print(f"  L1%p={L1%p} L2%p={L2%p} L3%p={L3%p}")
# quotient handles (only valid if L_i ≡0 mod p)
if L1%p==0: val[30317] = -(L1)//p
if L3%p==0: val[2936] = (537773*L3)//p
if L2%(6672769*p)==0: val[5146] = L2//(6672769*p)
else: print(f"  L2 %(6672769*p) = {L2%(6672769*p)} (NONzero -> x_5146 not integer-valid)")
H.forward()
F=recheck(); Fset=set(F)
print(f"\nAFTER residue-fix + quotient handles: sat={len(H.lines)-len(F)} nfail={len(F)}")
print(f"  core_fail={len(Fset & H.CORESET)} noncore_fail={len(Fset - H.CORESET)}")
print(f"  S%p={val[35389]%p} T%p={val[6671]%p}")
print(f"  first fails: {sorted(F)[:30]}")
if len(F)<r['nfail']:
    json.dump({f"x_{i}":val[i] for i in range(H.NVARS)}, open(f"agentD_quot_{quad}.json",'w'))
    print(f"  IMPROVED -> saved agentD_quot_{quad}.json ({len(H.lines)-len(F)})")
