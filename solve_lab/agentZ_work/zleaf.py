#!/usr/bin/env python3
"""Agent Z, step 1: does ANY pair of leaves share an x-coordinate?

Direct measurement from the 512 extracted leaf constants -- no assumption about
the leaf->exponent assignment.  Reduce each 287-296 bit constant mod p, pair
them per selector, orient (x,y) using the common curve relation, and look for
x-collisions.  A collision with equal y is the same point (a degeneracy);
a collision with unequal y is P's INFEASIBLE intermediate (N1 = -B^2 != 0).
"""
import os, json, collections, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'zsel.json')))
trip = D['triples']          # (selector, coord_var, K_string)
print("triples:", len(trip))

p = 2**256 - 2**32 - 977
Q = 97553848499418123410591666447050222001188385549510401465815187079080512838891
print("p == 2^256-2^32-977 :", p == 115792089237316195423570985008687907853269984665640564039457584007908834671663)

bysel = collections.defaultdict(list)
for s, cv, K in trip:
    bysel[s].append((cv, int(K) % p))
print("selectors with exactly 2 coordinates:", sum(1 for v in bysel.values() if len(v) == 2), "/", len(bysel))

# residues distinct?
res = [r for v in bysel.values() for cv, r in v]
print("distinct residues mod p among the 512 leaf constants:", len(set(res)))

# --- orientation: find the ordering making  y^2 - (x + Q/3)^3  a single common value
inv3 = pow(3, p - 2, p)
q3 = Q * inv3 % p

def bval(x, y):
    X = (x + q3) % p
    return (y * y - X * X % p * X) % p

sels = sorted(bysel)
cands = collections.Counter()
for s in sels:
    (c1, r1), (c2, r2) = bysel[s]
    cands[bval(r1, r2)] += 1          # order (x,y) = (r1,r2)
    cands[bval(r2, r1)] += 1
common, n = cands.most_common(1)[0]
print("most common curve constant b appears for", n, "of 512 orientations; b =", common)

leaves = {}
amb = 0
for s in sels:
    (c1, r1), (c2, r2) = bysel[s]
    o1 = bval(r1, r2) == common
    o2 = bval(r2, r1) == common
    if o1 and o2:
        amb += 1
    if o1:
        leaves[s] = (r1, r2)
    elif o2:
        leaves[s] = (r2, r1)
print("leaves oriented:", len(leaves), "/ 256   ambiguous:", amb)
assert len(leaves) == 256

# --- THE MEASUREMENT: x-collisions among the 256 leaves
xs = collections.defaultdict(list)
for s, (x, y) in leaves.items():
    xs[x].append(s)
dup = {x: v for x, v in xs.items() if len(v) > 1}
print()
print("distinct leaf x-coordinates:", len(xs), "/ 256")
print("LEAF PAIRS SHARING AN x-COORDINATE:", sum(len(v) * (len(v) - 1) // 2 for v in dup.values()))
for x, v in dup.items():
    ys = [leaves[s][1] for s in v]
    print("   x=%d  selectors=%s  ys equal:%s  ys negatives:%s"
          % (x, v, len(set(ys)) == 1, all((ys[0] + t) % p == 0 for t in ys[1:])))

# --- also: is any leaf the negative of another?  (P = -Q  <=>  same x)
negpairs = 0
xmap = {x: s for x, v in xs.items() for s in v}
for s, (x, y) in leaves.items():
    if x in xmap and xmap[x] != s:
        negpairs += 1
print("leaves that are +/- of another leaf:", negpairs)

json.dump({'p': str(p), 'q3': str(q3), 'b': str(common),
           'leaves': {str(s): [str(x), str(y)] for s, (x, y) in leaves.items()}},
          open(os.path.join(HERE, 'zleaves.json'), 'w'))
print("wrote zleaves.json")
