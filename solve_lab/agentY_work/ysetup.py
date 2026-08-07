#!/usr/bin/env python3
"""Agent Y -- build and verify the COMPLEMENT target T'.

Identity:   fold(S) + fold(Sbar) = sum_{i=0}^{255} 2^i G = (2^256 - 1) G  =: A
so if k*G = T with ON-set S, then k'*G = T' := A - T with ON-set complement(S),
and  w(k') = 256 - w(k).

Everything is re-derived from agent Q's raw instance-derived JSON (read-only) and
re-verified here; nothing is taken on trust.  Writes ydata.json.
"""
import json, os, sys, random

HERE = os.path.dirname(os.path.abspath(__file__))
Q    = os.path.join(HERE, '..', 'agentQ_work')

p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
cs = K * pow(3, p - 2, p) % p

cur = json.load(open(os.path.join(Q, 'curve.json')))
A_ = int(cur['a']); B_ = int(cur['b'])
assert int(cur['p']) == p and int(cur['c_shift']) == cs
print('p == 2^256-2^32-977 (secp256k1 prime) :', p == 2**256 - 2**32 - 977)
print('N == secp256k1 group order            :',
      N == 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141)
print('a, b                                  : a=%d  b=%d' % (A_, B_))

def inv(z): return pow(z, p - 2, p)
def add(P, Q_):
    if P is None: return Q_
    if Q_ is None: return P
    x1, y1 = P; x2, y2 = Q_
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0: return None
        l = (3 * x1 * x1 + A_) % p * inv(2 * y1 % p) % p
    else:
        l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
def neg(P): return None if P is None else (P[0], (-P[1]) % p)
def sub(P, Q_): return add(P, neg(Q_))
def mul(k, P):
    if k < 0: k, P = -k, neg(P)
    R = None; Qq = P
    while k > 0:
        if k & 1: R = add(R, Qq)
        Qq = add(Qq, Qq); k >>= 1
    return R
def oncur(P):
    return P is None or (P[1] * P[1] - pow(P[0], 3, p) - A_ * P[0] - B_) % p == 0

# ---------------------------------------------------------------- leaves / ladder
leaf = {int(g): (int(v[0]), int(v[1])) for g, v in
        json.load(open(os.path.join(Q, 'qleaf.json'))).items()}
print('leaves: %d, all on curve: %s' % (len(leaf), all(oncur(v) for v in leaf.values())))
lad = json.load(open(os.path.join(Q, 'qladder.json')))
e2s = {int(k): int(v) for k, v in lad['exp2sel'].items()}
G = leaf[e2s[0]]
cur_pt = G; ok = 0; bad = []
for i in range(256):
    if i in e2s:
        if leaf[e2s[i]] == cur_pt: ok += 1
        else: bad.append(i)
    cur_pt = add(cur_pt, cur_pt)
print('ladder L_i == 2^i*G                   : %d/%d  bad=%s' % (ok, len(e2s), bad))
assert not bad and ok == 256
assert mul(N, G) is None
print('N*G == O                              : True')

L = [G]
for i in range(255): L.append(add(L[-1], L[-1]))
assert all(oncur(P) for P in L)

# ---------------------------------------------------------------- target T
C1 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
C2 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
T = ((C1 + cs) % p, C2 % p)
assert oncur(T) and mul(N, T) is None
print('T on curve, N*T == O                  : True')
assert G == (31917591553801470078828036568057743875467637605644620066197178005619323650152,
             83364444556352143115103874010002344754157095926378075484791050960431190202517)
assert T == (30121525689829097248416773597728729849687459852468451992398421980273013515302,
             44859544763832475231923253825569092119321525945631045653619508440821028887)
print('G, T match the values Q and X searched: True')

# ---------------------------------------------------------------- A = (2^256-1) G, TWO WAYS
M = 2**256 - 1
# way 1: scalar multiplication (double-and-add) with the scalar reduced mod N
A1 = mul(M % N, G)
# way 2: fold every one of the 256 ladder leaves
A2 = None
for i in range(256):
    A2 = add(A2, L[i])
# way 3 (extra): scalar multiplication with the UNREDUCED 256-bit scalar
A3 = mul(M, G)
print()
print('A = (2^256-1)*G')
print('  way1 (scalar mult, M mod N)  = %s' % (A1,))
print('  way2 (fold all 256 leaves)   = %s' % (A2,))
print('  way3 (scalar mult, raw M)    = %s' % (A3,))
print('  way1 == way2                 : %s' % (A1 == A2))
print('  way1 == way3                 : %s' % (A1 == A3))
assert A1 == A2 == A3
A = A1
assert oncur(A) and mul(N, A) is None
print('  A on curve, N*A == O         : True')

# ---------------------------------------------------------------- T'
Tp = sub(A, T)
print()
print("T' = A - T")
print("  T'.x = %d" % Tp[0])
print("  T'.y = %d" % Tp[1])
print("  T' on curve                  : %s" % oncur(Tp))
print("  N*T' == O                    : %s" % (mul(N, Tp) is None))
print("  T + T' == A                  : %s" % (add(T, Tp) == A))
print("  T' != T                      : %s" % (Tp != T))
assert oncur(Tp) and mul(N, Tp) is None and add(T, Tp) == A

# ---------------------------------------------------------------- complement identity test
random.seed(20260807)
allbad = 0
for trial in range(12):
    w = random.randint(1, 40)
    S = sorted(random.sample(range(256), w))
    Sb = [i for i in range(256) if i not in set(S)]
    f1 = None
    for i in S: f1 = add(f1, L[i])
    f2 = None
    for i in Sb: f2 = add(f2, L[i])
    k = sum(1 << i for i in S)
    kb = sum(1 << i for i in Sb)
    good = (add(f1, f2) == A) and (k + kb == M) and (len(S) + len(Sb) == 256) \
           and (f1 == mul(k % N, G)) and (f2 == mul(kb % N, G))
    if not good: allbad += 1
print()
print('complement identity fold(S)+fold(Sbar)==A, k+kbar==2^256-1 : %d/12 trials ok'
      % (12 - allbad))
assert allbad == 0

# a target-level statement of exactly what the search means
print()
print("MEANING:  if k*G == T with ON-set S then k'*G == T' with ON-set complement(S),")
print("          and w(k') = 256 - w(k).  Exhausting w' <= W on T' with no hit proves")
print("          w' >= W+1, i.e. w <= 256 - (W+1) = 255 - W.")

json.dump({'p': str(p), 'a': str(A_), 'b': str(B_), 'N': str(N),
           'G': [str(G[0]), str(G[1])],
           'T': [str(T[0]), str(T[1])],
           'A': [str(A[0]), str(A[1])],
           'Tp': [str(Tp[0]), str(Tp[1])],
           'ladder': [[str(x), str(y)] for x, y in L]},
          open(os.path.join(HERE, 'ydata.json'), 'w'))
print('\nwrote ydata.json')
