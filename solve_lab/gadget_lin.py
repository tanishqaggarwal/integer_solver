#!/usr/bin/env python3
import heal_harness as H
import json
from math import gcd
p=H.p
d=json.load(open('fc_partial.json'))
fp={int(k[2:]) if k.startswith('x_') else int(k):int(v) for k,v in d.items()}
for v in H.freeinp: H.val[v]=fp.get(v,0)
H.forward()
F0=set(H.fails())
print(f"base: {len(F0)} fails")
# EXACT linear coeffs of x_31731 and x_9106 wrt x_14865, x_31861 (integer finite diff)
def exact_coef(target, w):
    b=H.val[target]; H.val[w]+=1; H.forward(); c=H.val[target]-b; H.val[w]-=1; H.forward(); return c
a1=exact_coef(31731,14865); a2=exact_coef(31731,31861)
b1=exact_coef(9106,14865)%13523997; b2=exact_coef(9106,31861)%13523997
x31731=H.val[31731]; x9106=H.val[9106]
print(f"x_31731={x31731}\n a1={a1}, a2={a2}")
# solve a1*d1 + a2*d2 = -x31731 (exact). g=gcd(a1,a2)
g=gcd(a1,a2)
print(f"gcd(a1,a2)={g}, x_31731 % g = {x31731 % g} (need 0 for exact solution)")
if x31731 % g == 0:
    # particular: use extended euclid
    def egcd(a,b):
        if b==0: return (a,1,0)
        gg,x,y=egcd(b,a%b); return (gg,y,x-(a//b)*y)
    gg,u,v=egcd(a1,a2)  # a1*u + a2*v = g
    k=(-x31731)//g
    d1_0=u*k; d2_0=v*k   # a1*d1_0 + a2*d2_0 = -x31731
    # homogeneous: d1 += t*(a2/g), d2 -= t*(a1/g)
    h1=a2//g; h2=-a1//g
    # now satisfy (B): (x9106 + b1*(d1_0+t*h1) + b2*(d2_0+t*h2)) % 13523997 == 0
    M=13523997
    base_b=(x9106 + b1*d1_0 + b2*d2_0)%M
    slope=(b1*h1 + b2*h2)%M
    print(f"mod-M: base={base_b}, slope={slope}, gcd(slope,M)={gcd(slope,M)}")
    if base_b % gcd(slope,M)==0 and slope%M!=0:
        t=(-base_b * pow(slope//gcd(slope,M), -1, M//gcd(slope,M))) % (M//gcd(slope,M))
        d1=d1_0 + t*h1; d2=d2_0 + t*h2
        print(f"solution: d1={d1 if abs(d1)<10**12 else 'big'}, d2={d2 if abs(d2)<10**12 else 'big'}")
        H.val[14865]+=d1; H.val[31861]+=d2; H.val[33168]=0
        H.forward()
        F1=set(H.fails())
        print(f"AFTER gadget solve: {len(F1)} fails; x_31731={H.val[31731]}, x_9106%13523997={H.val[9106]%13523997}")
        print(f"  fixed: {len(F0-F1)}, broke: {sorted(F1-F0)[:15]}")
        if len(F1)==0:
            json.dump({f'x_{i}':H.val[i] for i in range(H.NVARS)},open('sg_SOLVED.json','w')); print("*** SOLVED - saved sg_SOLVED.json ***")
    else:
        print("mod-M unsolvable with this slope")
else:
    print("x_31731 not divisible by gcd -> need more knobs")
