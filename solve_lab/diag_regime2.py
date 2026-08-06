#!/usr/bin/env python3
"""Check the nonlinear regime-2 core condition x_33469*x_1326^2 == x_27713^2 mod p (which decouples
x_29322,x_3558). Also verify the elimination is correct and check QR status."""
import json, re, sys
from propagate import NVARS
p=2**256-2**32-977
sol={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
def V(v): return sol.get(v,0)
x33469=V(33469); x1326=V(1326); x27713=V(27713)
x29322=V(29322); x3558=V(3558); x29356=V(29356); x35389=V(35389); x6671=V(6671)
print(f"x_33469 mod p = {x33469%p}")
print(f"x_1326   mod p = {x1326%p}")
print(f"x_27713  mod p = {x27713%p}")
# verify S = x_33469*x_29356 - x_3558^2  and  T = x_27713*x_29322 - x_3558*x_1326  (need x_11602,x_1326)
print(f"\nS=x_35389 mod p = {x35389%p}")
print(f"x_33469*x_29356 - x_3558^2 mod p = {(x33469*x29356 - x3558*x3558)%p}")
print(f"T=x_6671 mod p = {x6671%p}")
# regime-2 condition
lhs=(x33469*x1326*x1326)%p; rhs=(x27713*x27713)%p
print(f"\nREGIME-2: x_33469*x_1326^2 mod p = {lhs}")
print(f"          x_27713^2       mod p = {rhs}")
print(f"          gap = {(lhs-rhs)%p}  ({'HOLDS' if lhs==rhs else 'differs'})")
# QR check
def is_qr(a): return pow(a%p,(p-1)//2,p)==1
print(f"x_33469 is QR mod p: {is_qr(x33469)}")
print(f"(x_27713/x_1326)^2 mod p = {(x27713*pow(x1326,p-2,p))**2 %p if x1326%p else 'x1326=0'}")
# what would x_3558 be in regime 2: x_27713*x_29322/x_1326
if x1326%p:
    x3558_r2=(x27713*x29322*pow(x1326,p-2,p))%p
    print(f"regime-2 x_3558 target mod p = {x3558_r2}; current x_3558 mod p = {x3558%p}")
