# Agent J log

## t0 — baseline verified
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> satisfied 39026/39033, failing [12231, 12270, 12350, 14584, 18673, 22044, 29125]. CONFIRMED.

## Independent parse
Wrote jparse.py: uses Python `ast` on `x_N -> XN` rewrite; peels outer wrapper
(square / const multiplier / c1*S+c2*S) and decomposes S as left-nested
A0 + c1*A1 + c2*A2 ... chain.

## BREAKTHROUGH-CANDIDATE: the residual is 4 congruences in a 38-variable cone
Independent parse validated exactly (jvalidate: 0/39033 mismatches).
Model: each eq = mult * (sum c_i A_i)^k, k in {1,2,4}, mult != 0, 39033 distinct atoms,
all of degree <= 2.  Therefore the system is exactly `M a = 0` in the atom values.
At best/new_instance_partial_39026.json EXACTLY 7 atoms are nonzero (jatomval.py):
  a23326 (x7068 - x2099) - 7376877*x642
  a23327 x28730 - x17499*x9413        (x17499 = p)
  a35889 x29854 - x22665*x1329        (x22665 = p)
  a35890 5113045*x7075*x9118 - x29854
  a35891 x31864 - x28961*x10903       (x28961 = p)
  a35892 x7075*x8731 + x31864
  a35893 x642 - x28599*x17325         (x28599 = p)
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663 (secp256k1 p) CONFIRMED.
Residual == 4 conditions:
  (A) x4432  == x19964  (mod p)
  (B) x7068  == x2099   (mod 7376877*p)
  (C) p | (1-b1*b2)*x9118        (D) p | (1-b1*b2)*x8731
with b1=x2081, b2=x4287 free booleans and the 3-way MUX
  x2099  = b2(1-b1)*x31861 + b1(1-b2)*x6418  + b1b2*x9118
  x19964 = b2(1-b1)*x14865 + b1(1-b2)*x12553 + b1b2*x8731
Backward cone of the whole residual = 38 variables, 10 free inputs
  (x9118, x8731, x9413, x17325, x2081, x4287, x12553, x14865, x31861, x6418).
x_4432 and x_7068 have NO other definer -> they look like free inputs too.
NEXT: build a forward propagation engine and test setting x4432:=x19964, x7068:=x2099.

## Reduced mod-p model (jmodp.py) — the instance in 13 numbers
* Constraint atoms (non-definer) = 8458.  At the on-manifold base, exactly 4 are
  nonzero mod p: a8583, a30271, a35890, a35892.
* Only 26 constraints are moved by the 13 params + 2 booleans + x8731/x9118.
* Degrees measured exactly (finite differences over GF(p)): all chain constraints
  are DEGREE 1 in their param; only a20407/a20409/a31575 are degree 2-3.
* The chain is
    x6418<-a3895(pin C1)  x12553<-a3897(pin C2)  x22152<-a32257(pin C3)  x33462<-a32259(pin C4)
    x14853<-a30271  x24548<-a8583  x14623<-a22688  x31339<-a26603
    x8778<-a34370   x16742<-a27640 x22649<-a2694
    x22162<-a31571  x30213<-a731
  One Gauss-Seidel sweep solves it (jchain.py) -> only [20407,20409,31575] remain.
  CONFIRMS prior sessions' §129-131 independently.
* x15298 = x7715*x34554, x7715 = f(x2081), x34554 = f(x24601).  Setting EITHER
  boolean to 0 makes x15298 = 0, which kills a20407/a20409/a31575 outright.
* jsolve2.py (branch + targeted sweep) mod p:
    b1=0 b2=0 -> 3 violated [731, 24075, 31571]
    b1=0 b2=1 -> 2 violated [731, 31571]
    b1=1 b2=0 -> 2 violated [731, 31571]
    b1=1 b2=1 -> 3 violated [20407, 20409, 31575]
  a731  : x18956 == C5 (mod p)        a31571: x24468 == x13682 (mod p)
  Both are the OUTPUT pins, and with x15298 = 0 the coordinates cannot reach them.
