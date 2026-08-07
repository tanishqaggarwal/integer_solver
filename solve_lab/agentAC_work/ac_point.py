#!/usr/bin/env python3
"""Agent AC. Everything about T-as-a-point that is checkable from the instance
data alone and would, if it hit, be evidence against `popcount(log_G T)` being
distributed like the popcount of a uniform scalar.

NOT an investigation of how the instance was generated: every test below is a
question about a given point on a given curve, answerable without any reference
to a generator, a seed, or an emission order.
"""
import math, json, sys

# curve data: taken from agentX_work/xdata.json and agentY_work/ydata.json,
# which are INDEPENDENT parses of EQUATIONS.txt and agree bit-for-bit.
# Everything structural about them is re-verified here from scratch.
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
a = 0
b = 64019533680030876408443198762210829058751700634554282185987325820393598524794
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
Gx = 31917591553801470078828036568057743875467637605644620066197178005619323650152
Gy = 83364444556352143115103874010002344754157095926378075484791050960431190202517
Tx = 30121525689829097248416773597728729849687459852468451992398421980273013515302
Ty = 44859544763832475231923253825569092119321525945631045653619508440821028887

G = (Gx, Gy); T = (Tx, Ty)
O = None

def inv(x): return pow(x, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    if P == Q:
        l = (3 * x1 * x1 + a) * inv(2 * y1) % p
    else:
        l = (y2 - y1) * inv(x2 - x1) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
def mul(k, P):
    R = None; Q = P
    while k:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R
def oncurve(P):
    if P is None: return True
    x, y = P
    return (y * y - x * x * x - a * x - b) % p == 0

R = {}
def rec(name, val, note=""):
    R[name] = val
    print("%-52s %s   %s" % (name, val, note))

print("=" * 100)
print("0. STRUCTURE, re-verified here")
print("=" * 100)
import sympy
rec("p prime", sympy.isprime(p))
rec("N prime", sympy.isprime(N), "=> the ONLY subgroups are {O} and E(F_p); no small-index subgroup can exist")
rec("p == 2^256-2^32-977", p == 2**256 - 2**32 - 977)
rec("a == 0  (j-invariant 0, Aut = mu_6)", a == 0)
rec("G on curve", oncurve(G))
rec("T on curve", oncurve(T))
rec("N*G == O", mul(N, G) is None)
rec("N*T == O", mul(N, T) is None)
rec("T != O and T != G", T is not None and T != G)
# curve isomorphic to secp256k1 y^2=x^3+7 ?
u6 = b * inv(7) % p
rec("b/7 is a 6th power mod p", pow(u6, (p - 1) // 6, p) == 1,
    "=> F_p-isomorphic to secp256k1")

print()
print("=" * 100)
print("1. SUBGROUP / ORDER TESTS")
print("=" * 100)
rec("order(T) == N (forced: N prime, T != O)", sympy.isprime(N) and mul(N, T) is None and T is not None)
rec("T has no small-index subgroup to live in", True, "N prime: index of <T> is 1")
# x=0 points would be the order-3 points; they exist iff b is a QR
rec("b is a QR mod p (i.e. order-3 points exist)", pow(b, (p - 1) // 2, p) == 1,
    "must be False since 3 does not divide N")

print()
print("=" * 100)
print("2. AUTOMORPHISM TESTS (Aut(E) = mu_6 since j = 0)")
print("=" * 100)
# zeta_3 in F_p, beta with beta^3=1
zeta = pow(2, (p - 1) // 3, p)
if zeta == 1:
    for g in range(3, 100):
        zeta = pow(g, (p - 1) // 3, p)
        if zeta != 1: break
rec("zeta_3^3 == 1, zeta_3 != 1", pow(zeta, 3, p) == 1 and zeta != 1)
# the 6 automorphisms: (x,y) -> (zeta^i x, +-y)
auts = []
for i in range(3):
    z = pow(zeta, i, p)
    for s in (1, p - 1):
        auts.append((z, s))
fixed = []
for (z, s) in auts:
    if (z, s) == (1, 1): continue
    img = (z * Tx % p, s * Ty % p)
    if img == T: fixed.append((z, s))
rec("T fixed by some NON-TRIVIAL automorphism", len(fixed) > 0,
    "impossible unless x(T)=0 or y(T)=0")
rec("x(T) == 0", Tx == 0)
rec("y(T) == 0", Ty == 0)

# lambda = the scalar realised by (x,y)->(beta x, y).  Identify beta and lambda.
lam_cands = [r for r in sympy.nthroot_mod(1, 3, N, all_roots=True) if r != 1]
lam = None; beta = None
for z in [pow(zeta, 1, p), pow(zeta, 2, p)]:
    for L in lam_cands:
        if mul(L, G) == (z * Gx % p, Gy):
            lam, beta = L, z
rec("found (beta, lambda) with [lambda]P == (beta x, y)", lam is not None)
if lam is not None:
    ok = all(mul(lam, mul(2**i, G)) == (beta * mul(2**i, G)[0] % p, mul(2**i, G)[1]) for i in [0,1,5,17,100])
    rec("verified on 5 ladder points", ok)
    rec("lambda^3 == 1 mod N", pow(lam, 3, N) == 1)

print()
print("(section 3 -- the two 2^20 reachability sweeps -- is in ac_point2.py,\n which uses inversion-free Jacobian arithmetic)")

print()
print("=" * 100)
print("4. IS x(T) or y(T) A SPECIAL INTEGER?")
print("=" * 100)
for nm, v in [("x(T)", Tx), ("y(T)", Ty)]:
    rec(nm + " bit length", v.bit_length(), "(p is %d bits)" % p.bit_length())
    rec(nm + " popcount", bin(v).count("1"), "(uniform: 128 +- 8)")
    rec(nm + " NAF weight", bin((3 * v) ^ v).count("1"), "(uniform: ~85.3)")
    r = sympy.integer_nthroot(v, 2); rec(nm + " is a perfect square (over Z)", r[1])
    r = sympy.integer_nthroot(v, 3); rec(nm + " is a perfect cube (over Z)", r[1])
    for cn, c in [("0", 0), ("p", p), ("2^255", 2**255), ("2^256", 2**256)]:
        d = abs(v - c)
        rec("  |%s - %s| < 2^40" % (nm, cn), d < 2**40, "distance = 2^%.2f" % (math.log2(d) if d else 0))
    rec(nm + " < 2^40 (i.e. a small integer)", v < 2**40)
    rec(nm + " is a power of two", v != 0 and (v & (v - 1)) == 0)

print()
print("=" * 100)
print("5. MISCELLANEOUS CHEAP DISTINGUISHERS")
print("=" * 100)
rec("x(T) is a QR mod p", pow(Tx, (p - 1) // 2, p) == 1, "(50/50 under any model; recorded, not evidence)")
rec("y(T) is a QR mod p", pow(Ty, (p - 1) // 2, p) == 1, "(50/50)")
rec("y(T) is the even root (y < p/2)", Ty < p // 2)
rec("T == -G", T == (Gx, (-Gy) % p))
rec("T == G", T == G)
# is T a small multiple of ANY ladder point (equivalently 2^-i T small)?  covered by 3.
# does x(T) divide or relate to p trivially?
rec("gcd(x(T), p) != 1", math.gcd(Tx, p) != 1)
rec("x(T) == y(T)", Tx == Ty)
print()
print("=" * 100)
print("SUMMARY: number of tests that HIT (i.e. found structure):")
hits = [k for k, v in R.items() if v is True and k.startswith(("T in", "m*T", "T ==", "x(T) is a perfect",
        "y(T) is a perfect", "x(T) < ", "y(T) < ", "T fixed"))]
print("  structural hits:", hits if hits else "NONE")
json.dump({k: (v if isinstance(v, (bool, int, str)) else str(v)) for k, v in R.items()},
          open("ac_point.json", "w"), indent=1)
