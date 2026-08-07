#!/usr/bin/env python3
"""Agent Z, step 1b/2: recover the leaf->exponent map by DOUBLING (no assumption),
then settle P's infeasible-intermediate condition exactly.

Curve (P's law, re-derived in zleaf.py):  y^2 = X^3 + b  over F_p,  X = x + Q/3.
If the 256 leaves are {2^i * G}, then doubling each leaf lands on another leaf,
which is a decisive structural check costing 256 point doublings.
"""
import os, json, collections

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'zleaves.json')))
p = int(D['p']); q3 = int(D['q3']); b = int(D['b'])
leaves = {int(s): (int(x), int(y)) for s, (x, y) in D['leaves'].items()}
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
print("p bits", p.bit_length(), " N bits", N.bit_length())

# work in X = x + q3 coordinates
pts = {s: ((x + q3) % p, y) for s, (x, y) in leaves.items()}
for s, (X, y) in pts.items():
    assert (y * y - (X * X % p) * X - b) % p == 0, "leaf not on y^2=X^3+b"
print("all 256 leaves verified on  y^2 = X^3 + b  (mod p)")

def dbl(P):
    X, y = P
    lam = (3 * X * X % p) * pow(2 * y % p, p - 2, p) % p
    X2 = (lam * lam - 2 * X) % p
    y2 = (lam * (X - X2) - y) % p
    return (X2, y2)

byX = {}
for s, P in pts.items():
    byX.setdefault(P[0], []).append(s)
lookup = {P: s for s, P in pts.items()}

succ = {}
hits = 0
for s, P in pts.items():
    t = lookup.get(dbl(P))
    if t is not None:
        succ[s] = t
        hits += 1
print("leaves whose DOUBLE is another leaf: %d / 256" % hits)

# chain
pred = {t: s for s, t in succ.items()}
starts = [s for s in pts if s not in pred]
print("chain start(s) (not the double of any leaf):", starts)
if len(starts) == 1 and hits == 255:
    order = [starts[0]]
    while order[-1] in succ:
        order.append(succ[order[-1]])
    print("doubling chain length:", len(order))
    expo = {s: i for i, s in enumerate(order)}
    json.dump({str(s): i for s, i in expo.items()},
              open(os.path.join(HERE, 'zexpo.json'), 'w'))
    print("=> the 256 leaves are exactly {2^i * L0 : i = 0..255};  wrote zexpo.json")
else:
    expo = None
    print("=> leaves do NOT form a single doubling chain; exponent map not recovered")

# ---------------- step 1 (theory): 2^(i-j) == +-1 mod N ?
def factor(n):
    f = {}
    d = 2
    while d * d <= n and d < 10**7:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

known = [2, 3, 149, 631, 107361793816595537, 174723607534414371449,
         341948486974166000522343609283189]
m = N - 1
fac = {}
for q in known:
    while m % q == 0:
        fac[q] = fac.get(q, 0) + 1; m //= q
print("N-1 residual after known factors:", m)
if m != 1:
    fac.update(factor(m))
prod = 1
for q, e in fac.items():
    prod *= q ** e
print("factorisation of N-1 verified:", prod == N - 1, sorted(fac.items()))

order2 = N - 1
for q in fac:
    while order2 % q == 0 and pow(2, order2 // q, N) == 1:
        order2 //= q
print("ord_N(2) =", order2, " (bits %d)" % order2.bit_length())
print("ord_N(2) | (i-j) for some 1<=i-j<=255 ?", any(order2 <= 255 and d % order2 == 0 for d in range(1, 256)))
half = order2 // 2 if order2 % 2 == 0 else None
print("ord_N(2) even?", order2 % 2 == 0, " half =", half,
      " -> 2^d == -1 mod N needs d = half; half <= 255 ?", (half is not None and half <= 255))
print("DIRECT CHECK: any 1<=d<=255 with 2^d == +-1 mod N :",
      [d for d in range(1, 256) if pow(2, d, N) in (1, N - 1)])

# ---------------- step 2: the exact intermediate condition
print()
print("=" * 78)
print("INTERMEDIATE CONDITION")
print("A merge node v combines A = sum_{i in S n Lv} 2^i G and B = sum_{i in S n Rv} 2^i G,")
print("with Lv, Rv DISJOINT leaf supports.  x(A)=x(B) with y(A)!=y(B)  <=>  A = -B  <=>")
print("   k(S n Tv) := sum_{i in S n Tv} 2^i  ==  0  (mod N),   Tv = Lv u Rv.")
print("Because 0 <= k < 2^256 and 2^256 < 2N, the ONLY nonzero solution is k == N exactly.")
print("  2^256 < 2N :", 2**256 < 2 * N)
print("  popcount(N) =", bin(N).count('1'), " N bit positions span 0 ..", N.bit_length() - 1)
print("So a merge at v is infeasible for configuration S  IFF  S n Tv == bits(N) exactly,")
print("which requires bits(N) subset Tv, hence |Tv| >= popcount(N) =", bin(N).count('1'))
for size in (256, 178, 128, 78, 64):
    print("   a node with |Tv| = %3d can host it: %-5s ; if so it excludes 2^(256-|Tv|) = 2^%d subsets"
          % (size, size >= bin(N).count('1'), 256 - size))
