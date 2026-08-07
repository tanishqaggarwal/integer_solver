#!/usr/bin/env python3
"""The containment lattice, computed rather than argued.

Signed-digit distance  d(a,b) = w+-(a-b)  = min # of terms  +-2^e (distinct e, e<=255)
summing to (a-b) mod N, minimised over the two integer representatives in (-N,N).
This is a metric on Z/N (subadditive, symmetric, zero iff equal), and

    S(c,m) = { k : d(k,c) <= m }

is exactly the ball of radius m about c.  Two offsets' classes are disjoint at depth m
iff d(c_i,c_j) > 2m.  So the whole coverage question is a packing question.
"""
import json, os, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'aa_offsets.json')))
N = int(D['N'])
man = D['manifest']

def naf(n):
    out = 0; i = 0
    while n:
        if n & 1:
            z = 2 - (n & 3); out += 1; n -= z
        n >>= 1; i += 1
    return out
def nafmax(n):
    e = 0; i = 0
    while n:
        if n & 1:
            z = 2 - (n & 3); e = i; n -= z
        n >>= 1; i += 1
    return e
def w(x):
    """minimal signed weight with exponents <= 255, over both representatives"""
    r = x % N
    if r == 0: return 0
    best = 10**9
    for v in (r, r - N):
        a = abs(v)
        if nafmax(a) <= 255: best = min(best, naf(a))
    return best

M = 7   # depth actually run
tags = [x['tag'] for x in man]
cs = [int(x['c']) for x in man]
tier = {x['tag']: x['tier'] for x in man}

print('offsets: %d' % len(tags))
print('\n--- (L1) containment in the plain search:  S(c,m) subset S(0, m + reach(c)) ---')
print('%-11s %-6s %-6s %s' % ('tag', 'tier', 'reach', 'redundant at m<=7 ?'))
for x in sorted(man, key=lambda z: (z['tier'], z['reach'] if z['reach'] is not None else 9999)):
    r = x['reach']
    print('%-11s %-6d %-6s %s' % (x['tag'], x['tier'], r,
          'YES -- fully inside S(0,%d)' % (M + r) if r is not None and r <= 2 else
          ('partially (inside S(0,%d), itself unreachable)' % (M + r) if r is not None else 'n/a')))

print('\n--- (L5) pairwise signed-digit distances between offsets ---')
mind = 10**9; argmin = None
hist = {}
for i, j in itertools.combinations(range(len(cs)), 2):
    dd = w(cs[i] - cs[j])
    hist[dd] = hist.get(dd, 0) + 1
    if dd < mind: mind, argmin = dd, (tags[i], tags[j])
print('pairs                : %d' % (len(cs) * (len(cs) - 1) // 2))
print('minimum distance     : %d   between %s and %s' % (mind, argmin[0], argmin[1]))
print('pairs with d <= 2m=14: %d' % sum(v for k, v in hist.items() if k <= 2 * M))
print('smallest 12 distances:', sorted(hist.items())[:12])
print('\nDISJOINTNESS at m<=%d: %s' % (M,
      'ALL %d classes pairwise DISJOINT -- coverage is exactly %dx the plain search'
      % (len(cs), len(cs)) if mind > 2 * M else
      'NOT all disjoint; %d pairs overlap' % sum(v for k, v in hist.items() if k <= 2 * M)))

# size of one ball
from math import comb
tot = sum(comb(256, m) * 2**m for m in range(M + 1))
import math
print('\n|S(c,%d)| <= %d  ~ 2^%.1f    (union over %d offsets ~ 2^%.1f of the 2^256 keyspace)'
      % (M, tot, math.log2(tot), len(cs), math.log2(tot * len(cs))))
print('fraction of keyspace covered: 2^%.1f' % (math.log2(tot * len(cs)) - 256))

print('\n--- the overlapping pairs (d <= 14): which offsets are near-duplicates ---')
ov = []
for i, j in itertools.combinations(range(len(cs)), 2):
    dd = w(cs[i] - cs[j])
    if dd <= 2 * M: ov.append((dd, tags[i], tags[j]))
ov.sort()
for dd, a, b in ov: print('  d=%-3d %-11s %-11s' % (dd, a, b))
core = [i for i, x in enumerate(man) if x['tier'] <= 2 and x['tag'] not in
        ('2p256', 'n2p256', 'p2p256p1', 'n2p256m1', 'p2p257', 'n2p257', 'ones')]
mind2 = min(w(cs[i] - cs[j]) for i, j in itertools.combinations(core, 2))
print('\nEXCLUDING the tier-3 controls AND the exponent-256/257 cluster (%d offsets left):'
      '  min pairwise distance = %d  -> %s' % (len(core), mind2,
      'pairwise DISJOINT at m<=7' if mind2 > 2 * M else 'still overlapping'))

# connected components of the "overlap graph" = the honest number of distinct classes
par = list(range(len(cs)))
def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for dd, a, b in ov:
    ra, rb = find(tags.index(a)), find(tags.index(b))
    if ra != rb: par[ra] = rb
comp = {}
for i, t in enumerate(tags): comp.setdefault(find(i), []).append(t)
print('\n--- honest coverage multiplier: connected components of the d<=2m overlap graph ---')
print('%d offsets collapse to %d structurally distinct classes at m<=%d' % (len(cs), len(comp), M))
for k, v in sorted(comp.items(), key=lambda z: -len(z[1])):
    if len(v) > 1: print('   merged (%d): %s' % (len(v), ' '.join(v)))
import math
print('\nEFFECTIVE coverage = %d x |S(.,7)| ~ 2^%.1f  = 2^%.1f of the keyspace'
      % (len(comp), math.log2(tot * len(comp)), math.log2(tot * len(comp)) - 256))
