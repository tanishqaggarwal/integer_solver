#!/usr/bin/env python3
"""K28: the only remaining honest shots at a FULL solve.
The instance is now known to have a unique solution: the leaf ON-set must be the binary
support of the integer k with k*G = T in the chord-composition group.
 (a) does the group order N factor?  (small factors => the exponent splits)
 (b) is k small / near N / of special shape?  (bounded baby-step giant-step)
"""
import sys, os, json, time, random, math
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
import fold as FD
P = FD.P
N = int(json.load(open(K + '/order.json'))['N'][0])
D = FD.points()
ch = json.load(open(K + '/chain.json'))
bypow = {}
for i_s, e in ch['exp'].items():
    bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
G = bypow[0]
T = (int(D['target']['X']), int(D['target']['Y']))


def mr(n, k=40):
    if n < 2: return False
    for q in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n % q == 0: return n == q
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True


print('N =', N)
print('N prime (Miller-Rabin, 40 rounds):', mr(N))
sm = [q for q in range(2, 200000) if N % q == 0]
print('small prime factors < 2e5:', sm)


def mulk(k):
    R = FD.INF; i = 0
    while k:
        if k & 1: R = FD.add(R, bypow[i])
        k >>= 1; i += 1
    return R


# (b) bounded BSGS for k in [0, W) and in (N-W, N)
def bsgs(Q, base, W, label):
    m = int(math.isqrt(W)) + 1
    t0 = time.time()
    tbl = {}
    R = FD.INF
    for j in range(m):
        tbl[R] = j
        R = FD.add(R, base)
    mB = R                      # m*base
    negmB = (mB[0], (-mB[1]) % P)
    cur = Q
    for i in range(m + 1):
        if cur in tbl:
            k = i * m + tbl[cur]
            print('%s: FOUND k = %d  (%.0fs)' % (label, k, time.time() - t0))
            return k
        cur = FD.add(cur, negmB)
    print('%s: no k < %d  (%.0fs, m=%d)' % (label, W, time.time() - t0, m))
    return None


W = 1 << 40
k = bsgs(T, G, W, 'k in [0,2^40)')
if k is None:
    negT = (T[0], (-T[1]) % P)
    bsgs(negT, G, W, '-k in [0,2^40)  (i.e. k near N)')
