#!/usr/bin/env python3
"""Agent Z: verify by construction that the infeasible-intermediate condition
fires at exactly one configuration, S = bits(N), and nowhere else.

Uses the leaf points recovered in zleaf.py and the exponent map recovered by
doubling in zdouble.py -- no assumption imported from any other agent."""
import os, json, random

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'zleaves.json')))
p = int(D['p']); q3 = int(D['q3']); b = int(D['b'])
leaves = {int(s): (int(x), int(y)) for s, (x, y) in D['leaves'].items()}
expo = {int(s): i for s, i in json.load(open(os.path.join(HERE, 'zexpo.json'))).items()}
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
P0 = {expo[s]: ((leaves[s][0] + q3) % p, leaves[s][1]) for s in leaves}   # i -> 2^i * L0

def add(A, B):
    if A is None: return B
    if B is None: return A
    (x1, y1), (x2, y2) = A, B
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 % p) * pow(2 * y1 % p, p - 2, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)

def summ(bits):
    A = None
    for i in bits:
        A = add(A, P0[i])
    return A

bN = [i for i in range(256) if (N >> i) & 1]
print("popcount(N) =", len(bN), " min bit", min(bN), " max bit", max(bN))

# 1. the whole set sums to the identity  (k = N  =>  k*G = O)
print("sum over bits(N) of 2^i*L0  ==  identity :", summ(bN) is None)

# 2. ANY split of bits(N) into two nonempty disjoint parts gives A = -B
rng = random.Random(3)
allneg = True
for trial in range(12):
    rng.shuffle(bN)
    cut = rng.randrange(1, len(bN))
    A = summ(bN[:cut]); B = summ(bN[cut:])
    ok = (A is not None and B is not None and A[0] == B[0] and (A[1] + B[1]) % p == 0 and A[1] != B[1])
    allneg &= ok
print("12 random splits of bits(N): A and B share x and have opposite y (INFEASIBLE):", allneg)

# 3. control: random subsets of the SAME size never collide
hits = 0
for trial in range(200):
    S = rng.sample(range(256), len(bN))
    rng.shuffle(S)
    cut = rng.randrange(1, len(S))
    A = summ(S[:cut]); B = summ(S[cut:])
    if A is not None and B is not None and A[0] == B[0]:
        hits += 1
print("200 random weight-%d subsets, random split: x-collisions =" % len(bN), hits)

print()
print("CONCLUSION: k(S n Tv) == 0 (mod N) with 0 <= k < 2^256 < 2N forces k == N exactly.")
print("bits(N) has weight %d and spans positions 0..255, so a node can host it only if" % len(bN))
print("its leaf support Tv has |Tv| >= %d AND contains every bit of N." % len(bN))
print("For the measured root split 178 | 78 the ONLY such node is the root, where the")
print("condition reads S == bits(N) exactly: ONE configuration out of 2^256, density 2^-256.")
print("Worst case over ALL binary trees on 256 leaves: the nodes with |Tv| >= %d form a" % len(bN))
print("root path of at most %d nodes, so the excluded fraction is at most %d * 2^-%d < 2^-186."
      % (256 - len(bN) + 1, 256 - len(bN) + 1, len(bN)))
