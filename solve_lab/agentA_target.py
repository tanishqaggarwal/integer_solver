#!/usr/bin/env python3
"""Verify S,T algebra and compute the target x_24908 mod p that makes S=T=0
given baseline slack. Bits affect S,T ONLY through X=x_24908."""
import json
from agentA_harness import p, load_solution, forward

base = load_solution('best/new_instance_partial_39013.json'); forward(base)
g = lambda i: base[i]
gm = lambda i: base[i] % p

# verify exact gate algebra
print("=== verify exact algebra ===")
print("x_35389 == x_17702 - x_27762 :", g(35389) == g(17702) - g(27762))
print("x_17702 == x_33469*x_29356   :", g(17702) == g(33469)*g(29356))
print("x_27762 == x_3558^2          :", g(27762) == g(3558)**2)
print("x_29356 == x_29322^2         :", g(29356) == g(29322)**2)
print("x_3558  == x_24908 - x_16742 :", g(3558) == g(24908) - g(16742))
print("x_6671  == x_32680 - x_11602 :", g(6671) == g(32680) - g(11602))
print("x_32680 == x_27713*x_29322   :", g(32680) == g(27713)*g(29322))
print("x_11602 == x_3558*x_1326     :", g(11602) == g(3558)*g(1326))

# slack scalars (mod p)
a = gm(33469)   # x_33469
d = gm(29322)   # x_29322
u = gm(27713)   # x_27713
w = gm(1326)    # x_1326
k = gm(16742)   # x_16742
X0 = gm(24908)  # baseline x_24908
print(f"\nslack: a=x_33469 d=x_29322 u=x_27713 w=x_1326 k=x_16742 (all mod p)")

# S = a*d^2 - (X-k)^2 ; T = u*d - (X-k)*w   where X=x_24908
def Smod(X): return (a*d*d - (X-k)**2) % p
def Tmod(X): return (u*d - (X-k)*w) % p
print(f"\ncheck formulas at baseline X0:")
print(f"  Smod(X0)={Smod(X0)}  vs actual S mod p={gm(35389)}  match={Smod(X0)==gm(35389)}")
print(f"  Tmod(X0)={Tmod(X0)}  vs actual T mod p={gm(6671)}   match={Tmod(X0)==gm(6671)}")

# consistency condition (slack only): a*w^2 == u^2 mod p
lhs = (a*w*w) % p; rhs = (u*u) % p
print(f"\nconsistency a*w^2 ≡ u^2 mod p ? {lhs==rhs}   (lhs-rhs={((lhs-rhs))%p})")

# target from T: x_3558 = X-k ≡ u*d/w ; target X = k + u*d/w
winv = pow(w, p-2, p)
x3558_target = (u*d*winv) % p
Xtarget_T = (k + x3558_target) % p
print(f"\nfrom T: target x_3558 = u*d/w = {x3558_target}")
print(f"        target X=x_24908 = {Xtarget_T}")
print(f"        Smod(Xtarget_T)={Smod(Xtarget_T)}  Tmod(Xtarget_T)={Tmod(Xtarget_T)}")

# from S: x_3558 = ±d*sqrt(a)
def sqrt_modp(n):
    n%=p
    if n==0: return [0]
    if pow(n,(p-1)//2,p)!=1: return []
    # p ≡ 3 mod 4
    r=pow(n,(p+1)//4,p)
    return [r, p-r]
roots_a = sqrt_modp(a)
print(f"\nsqrt(a) exists: {len(roots_a)>0}")
for ra in roots_a:
    xt = (d*ra)%p
    print(f"  x_3558=±d*sqrt(a) = {xt}  -> X={ (k+xt)%p }  Smod={Smod((k+xt)%p)} Tmod={Tmod((k+xt)%p)}")

# Also: the CORE_REDUCTION target was x_3558=0, x_29322=0 (needs slack change). Our new target
# keeps slack at baseline and only moves X=x_24908.
json.dump({'a':a,'d':d,'u':u,'w':w,'k':k,'X0':X0,
           'x3558_target':x3558_target,'Xtarget':Xtarget_T,
           'consistency':(lhs==rhs)}, open('agentA_target.json','w'))
print("\nsaved agentA_target.json")
