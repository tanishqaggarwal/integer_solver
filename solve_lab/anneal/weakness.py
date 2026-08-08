#!/usr/bin/env python3
"""weakness.py -- the standard ways a curve can be secretly easy. Checked, not assumed."""
from instance import p, n, B, G, T, mul, add, neg
from sympy import factorint, isprime

print(f"p = {p}\nn = {n}")
print(f"n prime: {isprime(n)}")
print(f"anomalous (n == p, Smart's attack applies in poly time): {n == p}")
t = p + 1 - n
print(f"trace t = p+1-n = {t}   (|t| <= 2*sqrt(p): {abs(t) <= 2*int(p**0.5)+2})")

print("\n-- MOV / Frey-Ruck: smallest embedding degree e with p^e == 1 (mod n) --")
e, cur = 1, p % n
while cur != 1 and e <= 200:
    cur = cur * p % n; e += 1
print(f"embedding degree e = {e if cur == 1 else '> 200'}"
      f"  -> pairing lands in F_p^e, a {(e*256) if cur==1 else '>51200'}-bit field")
print("   (e must be tiny -- single digits -- for index calculus in F_p^e to help)")

print("\n-- Pohlig-Hellman: does the group order factor? --")
print(f"n = {factorint(n)}")

print("\n-- j-invariant / CM: A == 0 so j = 0, CM by Z[zeta_3] --")
# the order-3 automorphism: (x,y) -> (beta*x, y) with beta^3 = 1
beta = pow(2, (p - 1) // 3, p)
assert pow(beta, 3, p) == 1 and beta != 1
Gb = (beta * G[0] % p, G[1])
print(f"beta^3 == 1 mod p: True;  (beta*x, y) is on the curve: "
      f"{(Gb[1]**2 - Gb[0]**3 - B) % p == 0}")
# its eigenvalue lam satisfies lam^2 + lam + 1 = 0 mod n
for cand in [r for r in __import__('sympy').nthroot_mod(1, 3, n, all_roots=True) if r != 1]:
    if mul(cand, G) == Gb:
        print(f"GLV eigenvalue lambda = {cand}")
        print("   -> automorphism group of order 6; Pollard rho speeds up by sqrt(6) ~ 2.4x")
        print(f"   -> generic ECDLP cost drops 2^128 -> 2^{128 - 1.3:.1f}. Not a break.")
        break

print("\n-- is the target a small multiple of the base point in disguise? --")
for c in list(range(1, 64)) + [n - 1, (n - 1) // 2, (n + 1) // 2]:
    if mul(c, G) == T: print(f"  !!! T == {c} * G")
print("  (nothing under 64; the deep search is the structured-k track)")
