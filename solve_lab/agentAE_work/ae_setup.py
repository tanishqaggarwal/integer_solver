#!/usr/bin/env python3
"""agent AE -- independent re-verification of the curve data, and derivation of
every target / centre this thread uses.  Nothing downstream reads xdata.json.

Writes ae_data.json.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- curve
d = json.load(open(os.path.join(HERE, '..', 'agentX_work', 'xdata.json')))
p = int(d['p']); a = int(d['a']); b = int(d['b']); N = int(d['N'])
G = (int(d['G'][0]), int(d['G'][1]))
T = (int(d['T'][0]), int(d['T'][1]))
ladder = [(int(x), int(y)) for x, y in d['ladder']]

O = None
def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % p == 0: return None
        l = (3 * P[0] * P[0] + a) * inv(2 * P[1]) % p
    else:
        l = (Q[1] - P[1]) * inv(Q[0] - P[0]) % p
    x = (l * l - P[0] - Q[0]) % p
    return (x, (l * (P[0] - x) - P[1]) % p)
def mul(k, P):
    k %= N
    R = None; Q = P
    while k:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R
def oncurve(P):
    return P is None or (P[1] * P[1] - P[0] ** 3 - a * P[0] - b) % p == 0

checks = {}
checks['p == 2^256-2^32-977'] = (p == 2**256 - 2**32 - 977)
checks['N == secp256k1 order'] = (N == 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141)
checks['a == 0'] = (a == 0)
checks['G on curve'] = oncurve(G)
checks['T on curve'] = oncurve(T)
checks['N*G == O'] = (mul(N, G) is None)
checks['N*T == O'] = (mul(N, T) is None)
checks['ladder[0] == G'] = (ladder[0] == G)

# independent doubling chain -- ladder[i] must be 2^i G
P = G; bad = 0
for i in range(256):
    if ladder[i] != P: bad += 1
    P = add(P, P)
checks['ladder[i] == 2^i G for all i (independent doubling)'] = (bad == 0)
checks['all 256 ladder points on curve'] = all(oncurve(L) for L in ladder)

# cross-check G,T against agent Y's independently-derived file (different extraction path)
ypath = os.path.join(HERE, '..', 'agentY_work', 'ydata.json')
if os.path.exists(ypath):
    y = json.load(open(ypath))
    def gp(o, k):
        v = o[k]
        if isinstance(v, dict): return (int(v['x']), int(v['y']))
        return (int(v[0]), int(v[1]))
    try:
        checks['G matches agentY ydata.json'] = (gp(y, 'G') == G)
        checks['T matches agentY ydata.json'] = (gp(y, 'T') == T)
    except Exception as e:
        checks['agentY cross-check'] = 'SKIPPED: %s' % e

# ---------------------------------------------------------------- GLV lambda / beta
# beta = nontrivial cube root of 1 mod p, lam = nontrivial cube root of 1 mod N,
# paired so that (beta*x, y) == [lam](x,y).
def cube_roots(m):
    # m = 1 mod 3 for both p and N here
    assert m % 3 == 1
    out = []
    g = 2
    while len(out) < 2:
        r = pow(g, (m - 1) // 3, m)
        if r != 1 and r not in out: out.append(r)
        g += 1
    return out
betas = cube_roots(p); lams = cube_roots(N)
beta = lam = None
for bb in betas:
    for ll in lams:
        ok = True
        for i in (1, 5, 37, 200):
            P = ladder[i]
            if ((bb * P[0]) % p, P[1]) != mul(ll, P): ok = False; break
        if ok: beta, lam = bb, ll
checks['GLV (beta,lam) found and verified on 4 ladder points'] = (beta is not None)
checks['lam^3 == 1 mod N'] = (pow(lam, 3, N) == 1 and lam != 1)
checks['beta^3 == 1 mod p'] = (pow(beta, 3, p) == 1 and beta != 1)
# verify phi == [lam] on 8 further random-ish points
r = 12345678901234567890
bad2 = 0
for i in range(8):
    P = mul(pow(r, i + 1, N) % N, G)
    if ((beta * P[0]) % p, P[1]) != mul(lam, P): bad2 += 1
checks['phi == [lam] on 8 more points'] = (bad2 == 0)

# ---------------------------------------------------------------- report
allok = all(v is True for v in checks.values())
for k, v in checks.items():
    print('%-55s %s' % (k, v))
print('ALL CHECKS PASS:', allok)
if not allok:
    sys.exit(1)

out = dict(p=str(p), a=str(a), b=str(b), N=str(N),
           G=[str(G[0]), str(G[1])], T=[str(T[0]), str(T[1])],
           beta=str(beta), lam=str(lam),
           ladder=[[str(x), str(y)] for x, y in ladder])
json.dump(out, open(os.path.join(HERE, 'ae_data.json'), 'w'))
print('wrote ae_data.json')
