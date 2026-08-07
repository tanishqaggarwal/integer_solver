#!/usr/bin/env python3
"""agent AF, step 12: the gate algebra.  L_v = OR(I_v) AND OR(J_v) ?  laminarity?  live count?"""
import sys, os, pickle, random
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
G = pickle.load(open(os.path.join(HERE, 'af_gates.pkl'), 'rb'))
gf = G['gatefn']; sel = G['sel']

def orleaves(n):
    """if n is a pure OR-tree over distinct leaves, return the frozenset of leaf indices;
       else return None"""
    if n[0] == 'leaf':
        return frozenset([n[1]])
    if n[0] == 'or':
        a = orleaves(n[1]); b = orleaves(n[2])
        if a is None or b is None:
            return None
        return a | b
    return None

pure = {}
bad = []
for L, f in gf.items():
    if f[0] == 'const':
        continue
    assert f[0] == 'and', f[0]
    I, J = orleaves(f[1]), orleaves(f[2])
    if I is None or J is None:
        bad.append(L); continue
    pure[L] = (I, J)
print('AND gates: %d ;  of the form OR(I) AND OR(J) with I,J pure OR-trees: %d ;  other: %d'
      % (sum(1 for f in gf.values() if f[0] == 'and'), len(pure), len(bad)))
dis = sum(1 for I, J in pure.values() if I & J)
print('gates with I ∩ J ≠ ∅ : %d' % dis)
print('|I|+|J| histogram (top):', sorted(Counter(len(I) + len(J) for I, J in pure.values()).items())[:8], '... max',
      max(len(I) + len(J) for I, J in pure.values()))

# ---- laminarity of the family {I_v} ∪ {J_v} ----
fam = set()
for I, J in pure.values():
    fam.add(I); fam.add(J)
print('distinct slot supports: %d' % len(fam))
faml = sorted(fam, key=len)
viol = 0
for i in range(len(faml)):
    for j in range(i + 1, len(faml)):
        A, B = faml[i], faml[j]
        if A & B and not (A <= B or B <= A):
            viol += 1
print('laminarity violations: %d' % viol)
roots = [x for x in fam if not any(x < y for y in fam)]
print('maximal sets: %d  sizes %s' % (len(roots), sorted(len(r) for r in roots)))
allleaf = set()
for r in roots:
    allleaf |= r
print('union of maximal sets covers %d of 256 selectors' % len(allleaf))
sz = Counter(len(x) for x in fam)
print('slot-support size histogram: %s' % sorted(sz.items())[:10])

# the parent of each slot support: smallest strict superset
byset = {}
for L, (I, J) in pure.values() and pure.items():
    byset[(I, J)] = L

# ---- live-count law ----
def live(sset):
    n = 0
    for L, (I, J) in pure.items():
        if (I & sset) and (J & sset):
            n += 1
    return n

print()
random.seed(1)
rows = []
for m in [0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 128, 144, 200, 233, 255, 256]:
    for t in range(3):
        S = frozenset(random.sample(range(256), m))
        rows.append((m, live(S), live(S) - max(m - 1, 0)))
print(' |S|   #live blocks   #live-(|S|-1)')
for m, l, d in rows:
    print('%4d %10d %12d' % (m, l, d))
pickle.dump({'pure': pure, 'fam': fam}, open(os.path.join(HERE, 'af_lam.pkl'), 'wb'), 2)
