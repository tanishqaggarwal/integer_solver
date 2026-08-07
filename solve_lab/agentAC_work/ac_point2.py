#!/usr/bin/env python3
"""Agent AC, part 5 continued: the O(2^20) reachability sweeps, done in Jacobian
coordinates so no modular inversion is needed in the inner loop.

Question asked (legitimate: it is a question about a given point):
  is T equal to  +- zeta_3^j * (c G)  for any tiny c?   [= Aut(E) o tiny scalar]
  is T equal to  +- zeta_3^j * (G / m) for any tiny m?  [i.e. m T = +- G]
Both are 'no' unless printed otherwise.
"""
import math, json, sys

p  = 115792089237316195423570985008687907853269984665640564039457584007908834671663
b  = 64019533680030876408443198762210829058751700634554282185987325820393598524794
N  = 115792089237316195423570985008687907852837564279074904382605163141518161494337
Gx = 31917591553801470078828036568057743875467637605644620066197178005619323650152
Gy = 83364444556352143115103874010002344754157095926378075484791050960431190202517
Tx = 30121525689829097248416773597728729849687459852468451992398421980273013515302
Ty = 44859544763832475231923253825569092119321525945631045653619508440821028887

LIMIT = 1 << 20

def jdouble(P):
    X1, Y1, Z1 = P
    if Z1 == 0 or Y1 == 0: return (0, 1, 0)
    A = X1*X1 % p; B = Y1*Y1 % p; C = B*B % p
    D = 2*((X1+B)*(X1+B) - A - C) % p
    E = 3*A % p; F = E*E % p
    X3 = (F - 2*D) % p
    Y3 = (E*(D - X3) - 8*C) % p
    Z3 = 2*Y1*Z1 % p
    return (X3, Y3, Z3)

def jadd_affine(P, q):
    X1, Y1, Z1 = P
    x2, y2 = q
    if Z1 == 0: return (x2, y2, 1)
    Z1Z1 = Z1*Z1 % p
    U2 = x2*Z1Z1 % p
    S2 = y2*Z1 % p * Z1Z1 % p
    if U2 == X1:
        if S2 == Y1: return jdouble(P)
        return (0, 1, 0)
    H = (U2 - X1) % p; HH = H*H % p; I = 4*HH % p; J = H*I % p
    r = 2*(S2 - Y1) % p; V = X1*I % p
    X3 = (r*r - J - 2*V) % p
    Y3 = (r*(V - X3) - 2*Y1*J) % p
    Z3 = ((Z1+H)*(Z1+H) - Z1Z1 - HH) % p
    return (X3, Y3, Z3)

def jx(P):
    X, Y, Z = P
    if Z == 0: return None
    zi = pow(Z, p-3, p)          # Z^-2 = Z^(p-3) since Z^(p-1)=1
    return X * zi % p

# --- sanity: Jacobian engine against a naive affine engine -------------------
def inv(x): return pow(x, p-2, p)
def aadd(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1 = P; x2,y2 = Q
    if x1 == x2 and (y1+y2) % p == 0: return None
    l = ((3*x1*x1)*inv(2*y1) if P == Q else (y2-y1)*inv(x2-x1)) % p
    x3 = (l*l - x1 - x2) % p
    return (x3, (l*(x1-x3) - y1) % p)
A = None; J = (0,1,0)
for i in range(1, 300):
    A = aadd(A, (Gx, Gy)); J = jadd_affine(J, (Gx, Gy))
    assert jx(J) == A[0], ("jacobian engine mismatch at", i)
print("jacobian engine == affine engine for c = 1..299   OK")

# zeta_3 in F_p
zeta = pow(2, (p-1)//3, p)
assert pow(zeta, 3, p) == 1 and zeta != 1
xs_T = {Tx % p, zeta*Tx % p, zeta*zeta % p * Tx % p}
xs_G = {Gx % p, zeta*Gx % p, zeta*zeta % p * Gx % p}
print("target x-set |{zeta^j x(T)}| =", len(xs_T), "   |{zeta^j x(G)}| =", len(xs_G))
print("(x-coordinate equality covers BOTH signs, since x(-P) = x(P))")

def sweep(base, targets, label):
    P = (0, 1, 0)
    hit = None
    for c in range(1, LIMIT+1):
        P = jadd_affine(P, base)
        X, Y, Z = P
        if Z == 0: continue
        ZZ = Z*Z % p
        for t in targets:
            if X == t*ZZ % p:
                hit = c; break
        if hit: break
        if c % 200000 == 0:
            print("   ... %s c=%d" % (label, c)); sys.stdout.flush()
    return hit

h1 = sweep((Gx, Gy), xs_T, "cG vs T")
print("T in { +- zeta^j (c G) : c <= 2^20 } :", "HIT c=%d" % h1 if h1 else "NO")
h2 = sweep((Tx, Ty), xs_G, "mT vs G")
print("G in { +- zeta^j (m T) : m <= 2^20 } :", "HIT m=%d" % h2 if h2 else "NO")

# the 1536-point endomorphism orbit  +- lambda^j 2^i G  (lambda acts as x -> beta x)
orb = False
P = (Gx, Gy)
Pj = (Gx, Gy, 1)
for i in range(256):
    X, Y, Z = Pj
    ZZ = Z*Z % p
    for t in xs_T:
        if X == t*ZZ % p: orb = True
    Pj = jdouble(Pj)
print("T in the 1536-point orbit { +- lambda^j 2^i G } :", orb)

json.dump({"limit": LIMIT, "hit_cG": h1, "hit_mT": h2, "orbit": orb},
          open("ac_point2.json", "w"))
