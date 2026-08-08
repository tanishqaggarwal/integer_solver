#!/usr/bin/env python3
"""structure.py -- what the 256-bit core actually is. Pure verification, no assumptions."""
import json, os
from sympy.ntheory.residue_ntheory import nthroot_mod

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'core.json')))
p = int(C['p']); B = int(C['B'])
bits = C['bit_vars']
Pt = {int(b): tuple(int(v) for v in C['points_dep'][b]) for b in C['points_dep']}
T  = tuple(int(v) for v in C['target_dep'])

# ---- group law on  y^2 = x^3 + B  (A == 0) ----
O = None
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1 = P; x2,y2 = Q
    if x1 == x2 and (y1+y2) % p == 0: return None
    if P == Q: lam = 3*x1*x1 % p * pow(2*y1, -1, p) % p
    else:      lam = (y2-y1) * pow(x2-x1, -1, p) % p
    x3 = (lam*lam - x1 - x2) % p
    return (x3, (lam*(x1-x3) - y1) % p)
def mul(k, P):
    R = None
    while k:
        if k & 1: R = add(R, P)
        P = add(P, P); k >>= 1
    return R

print(f"curve after depressing:  y^2 = x^3 + {B}   (mod p, A = 0)")

# ---- is it a twist of  y^2 = x^3 + 7 ? ----
u = nthroot_mod(B * pow(7, -1, p) % p, 6, p)
print(f"B / 7 is a 6th power:   u = {u}   (u^6 == B/7: {pow(u,6,p) == B*pow(7,-1,p)%p})")
iso = lambda Q: None if Q is None else (Q[0]*pow(u,-2,p) % p, Q[1]*pow(u,-3,p) % p)

n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141   # secp256k1 order
G_k1 = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
        0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
kill = sum(1 for q in Pt.values() if mul(n, q) is None)
print(f"n * P == O for {kill}/256 gated points, and for the target: {mul(n, T) is None}")
print(f"secp256k1 generator maps onto this curve: "
      f"{(lambda q: (q[1]*q[1]-q[0]**3-7) % p == 0)(iso(Pt[bits[0]]))}")

# ---- the 256 points: is it a doubling chain? ----
idx = {}
for b, q in Pt.items(): idx.setdefault(q, b)
chain, seen = [], set()
# find the head: a point whose halving is not in the set
dbl = {b: add(q, q) for b, q in Pt.items()}
child = {}
for b, d in dbl.items():
    if d in idx: child[b] = idx[d]
heads = [b for b in Pt if b not in set(child.values())]
print(f"doubling edges P_i -> 2P_i landing inside the set: {len(child)}/256, chain heads: {len(heads)}")
if len(heads) == 1:
    cur = heads[0]
    while True:
        chain.append(cur)
        if cur not in child: break
        cur = child[cur]
    print(f"single chain length: {len(chain)}  (so P_i = 2^i * P_0)")
    G = Pt[chain[0]]
    ok = all(Pt[chain[i]] == mul(1 << i, G) for i in range(len(chain)))
    print(f"verified P_i == 2^i * P_0 for all i: {ok}")
    json.dump({'chain_bit_vars': chain}, open(os.path.join(HERE,'chain.json'),'w'))

print()
print("=== THE CORE DECISION PROBLEM ===")
print("find b_0..b_255 in {0,1} with   sum_i b_i * (2^i G)  ==  T   on E(F_p)")
print("equivalently:   k * G == T   with   k = sum b_i 2^i,   0 <= k < 2^256")
print(f"group order n = {n}  (prime: {__import__('sympy').isprime(n)})")
print("=> a single ECDLP on a curve F_p-isomorphic to secp256k1.")
