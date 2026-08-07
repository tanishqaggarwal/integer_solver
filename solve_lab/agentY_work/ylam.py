#!/usr/bin/env python3
"""Agent Y -- the order-6 automorphism orbit, and its complements.

secp256k1 has phi(x,y) = (beta*x, y) with beta^3 = 1 mod p, and phi = multiplication by
lambda, lambda^3 = 1 mod N.  If k*G = T then

    phi(T)   = (lambda   * k) G
    phi^2(T) = (lambda^2 * k) G
    -T       = (N - k) G

so the SIX targets  {+-T, +-phi(T), +-phi^2(T)}  correspond to the six scalars
{+-k, +-lambda k, +-lambda^2 k} -- different integers, same difficulty class.  A low-weight
hit on ANY of them recovers k by a known automorphism.  Each also has a complement target
A - X (A = (2^256-1)G), giving TWELVE targets in total.

This script derives beta and lambda, verifies phi == [lambda] on random points, and writes
one engine input file per target.
"""
import json, os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'ydata.json')))
p = int(d['p']); A_ = int(d['a']); N = int(d['N'])
L = [(int(x), int(y)) for x, y in d['ladder']]
G = (int(d['G'][0]), int(d['G'][1]))
T = (int(d['T'][0]), int(d['T'][1]))
A = (int(d['A'][0]), int(d['A'][1]))

def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0: return None
        l = (3 * x1 * x1 + A_) % p * inv(2 * y1 % p) % p
    else:
        l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
def neg(P): return None if P is None else (P[0], (-P[1]) % p)
def sub(P, Q): return add(P, neg(Q))
def mul(k, P):
    k %= N
    R = None; Q = P
    while k > 0:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R

# --- beta: nontrivial cube root of 1 mod p
assert (p - 1) % 3 == 0
beta = None
g = 2
while beta is None:
    c = pow(g, (p - 1) // 3, p)
    if c != 1: beta = c
    g += 1
assert pow(beta, 3, p) == 1 and beta != 1
# --- lambda: nontrivial cube root of 1 mod N
assert (N - 1) % 3 == 0
lam = None
g = 2
while lam is None:
    c = pow(g, (N - 1) // 3, N)
    if c != 1: lam = c
    g += 1
assert pow(lam, 3, N) == 1 and lam != 1

def phi(P): return None if P is None else ((beta * P[0]) % p, P[1])

# which of the two lambdas matches this beta?
if phi(G) != mul(lam, G):
    lam = pow(lam, 2, N)
assert phi(G) == mul(lam, G)
print('beta   = %d' % beta)
print('lambda = %d' % lam)
print('beta^3 == 1 mod p : True     lambda^3 == 1 mod N : True')
random.seed(11)
ok = 0
for _ in range(8):
    r = random.randrange(1, N)
    P = mul(r, G)
    if phi(P) == mul(lam, P): ok += 1
print('phi(P) == [lambda] P on random points : %d/8' % ok)
assert ok == 8

lam2 = pow(lam, 2, N)
targets = {
    'T':        T,
    'negT':     neg(T),
    'lamT':     phi(T),
    'neglamT':  neg(phi(T)),
    'lam2T':    phi(phi(T)),
    'neglam2T': neg(phi(phi(T))),
}
print('\ntarget            scalar it encodes             on curve  N*X==O')
scal = {'T': 'k', 'negT': '-k', 'lamT': 'lam*k', 'neglamT': '-lam*k',
        'lam2T': 'lam^2*k', 'neglam2T': '-lam^2*k'}
allt = {}
for nm, P in targets.items():
    oc = (P[1] * P[1] - pow(P[0], 3, p) - A_ * P[0] - int(d['b'])) % p == 0
    print('  %-9s %-28s %-9s %s' % (nm, scal[nm], oc, mul(N, P) is None))
    allt[nm] = P
    allt['c_' + nm] = sub(A, P)        # complement:  (2^256-1) - that scalar
# distinctness
vals = list(allt.values())
print('\n12 targets, all distinct : %s' % (len(set(vals)) == 12))

lad = ''.join('%d %d\n' % (x, y) for x, y in L)
for nm, P in allt.items():
    with open(os.path.join(HERE, 'data_%s.txt' % nm), 'w') as f:
        f.write('%d %d\n' % (P[0], P[1])); f.write(lad)
print('wrote 12 engine input files data_<name>.txt')
json.dump({nm: [str(P[0]), str(P[1])] for nm, P in allt.items()},
          open(os.path.join(HERE, 'yorbit.json'), 'w'))
