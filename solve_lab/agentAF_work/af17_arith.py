#!/usr/bin/env python3
"""agent AF, step 17: the arithmetic content of the 927.  Does c divide the row coefficients?"""
import sys, os, pickle
from math import gcd
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
from af1_parse import is_const
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
conds = C['conds']; info = M['info']

def term(t):
    if t[0] == '*':
        k = is_const(t[1])
        if k is not None:
            return (k, t[2])
        k = is_const(t[2])
        if k is not None:
            return (k, t[1])
    if t[0] == 'neg':
        k, x = term(t[1]); return (-k, x)
    return (1, t)

rows = []
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    if info[i]['cls'] != 'cong':
        continue
    rhs = defc[info[i]['other']][0][1]
    (a1, D1), (a2, D2) = term(rhs[1]), term(rhs[2])
    rows.append((info[i]['gate'], c, a1, a2, D1, D2))
print('congruence rows: %d' % len(rows))
g1 = Counter()
for (g, c, a1, a2, D1, D2) in rows:
    if c == 1:
        continue
    g1[(gcd(c, a1) == c, gcd(c, a2) == c)] += 1
print('c>1 rows:  (c|alpha, c|beta) census:', dict(g1))
print('  gcd(c,alpha) distribution (c>1):',
      Counter(gcd(c, a1) for (g, c, a1, a2, D1, D2) in rows if c > 1).most_common(5))
print('  gcd(c, gcd(alpha,beta)) == c  for:',
      sum(1 for (g, c, a1, a2, D1, D2) in rows if c > 1 and gcd(c, gcd(a1, a2)) == c),
      'of', sum(1 for r in rows if r[1] > 1))

# factorisation of the c's
def fac(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

cs = [c for (R, c, Ex, aid, uc) in conds if c > 1]
print('\n927 multipliers: min %d max %d ; all distinct: %s' % (min(cs), max(cs), len(set(cs)) == len(cs)))
sq = 0; pr = Counter(); nf = Counter()
allp = Counter()
for c in cs:
    f = fac(c)
    nf[sum(f.values())] += 1
    for p, e in f.items():
        allp[p] += 1
    if any(e > 1 for e in f.values()):
        sq += 1
print('  #prime factors (with multiplicity) histogram:', dict(sorted(nf.items())))
print('  non-squarefree: %d' % sq)
print('  distinct primes across all 927: %d ; most common:' % len(allp), allp.most_common(8))
sm = sorted(p for p in allp)
print('  smallest primes used:', sm[:12])

# per-block: are the three c_k coprime?
byblk = defaultdict(list)
for (g, c, a1, a2, D1, D2) in rows:
    byblk[g].append(c)
pc = Counter()
for g, l in byblk.items():
    l = [c for c in l if c > 1]
    pc[len(l)] += 1
print('\nper block, #congruence rows with c>1:', dict(sorted(pc.items())))
pickle.dump({'rows': rows}, open(os.path.join(HERE, 'af_rows.pkl'), 'wb'), 2)
