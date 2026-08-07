#!/usr/bin/env python3
"""Planted-answer validation, per offset class, EXACT and end-to-end.

For each plant we know k = c + delta with delta a 5-term signed rep.  We predict, in pure
Python with an independent group law, the precise `HIT <sz> <code> <s_last> <key>` line the
engine must print for the split |alpha|=3 / |beta|=2 (and also |alpha|=2 / |beta|=3), then
check the engine printed it, then DECODE that line back to a scalar and verify k*G == T'
on the curve.

Algebra (this is where the two new failure modes live):
  scan side  Q = B + sum_{s in beta} sigma_s 2^{e_s} G,   B = (k - c) G = delta*G
  table side x(V*G),  V = sum_{alpha} eps 2^e, lowest-exponent eps forced to +1
  a match is x(Q) == x(V*G)  <=>  Q = +-V*G  <=>  delta + sum_beta = +- V.
Because only x is compared, the global-negation WLOG on the table is absorbed by the +-
branch: the scan indices are ALWAYS the negation of the complementary delta terms,
independent of which sign the table's leading term ended up with.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'aa_offsets.json')))
p = int(D['p']); N = int(D['N']); A_ = int(D['a'])
G = (int(D['G'][0]), int(D['G'][1]))
X = os.path.join(HERE, '..', 'agentX_work')
lad = [(int(a), int(b)) for a, b in json.load(open(os.path.join(X, 'xdata.json')))['ladder']]

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
def mul(k, P):
    k %= N; R = None; Q = P
    while k > 0:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R
def sp(s):
    P = lad[s >> 1]
    return P if not (s & 1) else (P[0], (-P[1]) % p)
def idx(e, sgn): return 2 * e + (0 if sgn > 0 else 1)
M64 = (1 << 64) - 1

def hits_of(path):
    out = set()
    if not os.path.exists(path): return out
    for ln in open(path):
        f = ln.split()
        if f and f[0] == 'HIT': out.add((int(f[1]), int(f[2]), int(f[3]), int(f[4])))
    return out

allok = True
print('%-8s %-6s %-6s %-9s %-9s %s' % ('plant', 'hits', 'split', 'lineseen', 'decodeOK', 'k unsigned wt'))
for pl in D['plants']:
    tag = pl['tag']; c = int(pl['c']); k = int(pl['k']); delta = int(pl['delta'])
    terms = sorted(zip(pl['exps'], pl['sgn']))
    base = (int(pl['base'][0]), int(pl['base'][1]))
    H = hits_of(os.path.join(HERE, 'runs', 'r_p_%s.txt' % tag))
    Tp = mul(k, G)
    assert base == mul(delta % N, G), 'base != (k-c)G for %s' % tag

    for na in (3, 2):                       # |alpha| = 3 and |alpha| = 2 splits
        alpha, beta = terms[:na], terms[na:]
        V0 = sum(s << e for e, s in alpha)
        g = 1 if alpha[0][1] > 0 else -1    # table stores V = g*V0 (lowest-exp sign forced +)
        V = g * V0
        bidx = sorted(idx(e, -s) for e, s in beta)          # negated complementary terms
        Q = base
        for s in bidx: Q = add(Q, sp(s))
        sz = len(bidx)
        code = 0
        for i, s in enumerate(bidx): code |= s << (16 * i)
        line = (sz, code, bidx[-1], Q[0] & M64)
        seen = line in H

        # decode the line back to a scalar, independently
        b = [(code >> (16 * i)) & 0xFFFF for i in range(sz)]
        Sb = sum((-1 if (s & 1) else 1) << (s >> 1) for s in b)
        dec = None
        for br in (1, -1):
            cand = (br * V - Sb) % N
            if mul(cand + c, G) == Tp:
                dec = cand; break
        okd = (dec is not None and (dec - delta) % N == 0
               and (mul(V % N, G)[0] & M64) == line[3])
        print('%-8s %-6d |a|=%-3d %-9s %-9s %d' %
              (tag, len(H), na, seen, okd, bin(k).count('1')))
        if not (seen and okd): allok = False

print('\nPLANT VALIDATION:',
      'PASS -- every offset class recovers its planted m=5 answer at both splits, '
      'and every recovered line decodes to the planted k with k*G == T on the curve'
      if allok else 'FAIL')
sys.exit(0 if allok else 1)
