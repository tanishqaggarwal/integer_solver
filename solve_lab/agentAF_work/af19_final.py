#!/usr/bin/env python3
"""agent AF, step 19: closing verifications.
   (1) the 180 X==Y diff conditions really are alias-identical;
   (2) the 3x2 congruence coefficient matrix has rank 2 for every block;
   (3) the gate family really is a binary tree -> the live-count law is a theorem;
   (4) the exact |S|-profile of every condition family.
"""
import sys, os, pickle, itertools, random
from math import gcd
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
R = pickle.load(open(os.path.join(HERE, 'af_rows.pkl'), 'rb'))
LAM = pickle.load(open(os.path.join(HERE, 'af_lam.pkl'), 'rb'))
conds = C['conds']; info = M['info']; rows = R['rows']; pure = LAM['pure']

# ---- (1) alias-identical diffs ----
copyset = set()
for aid, a in enumerate(atoms):
    if a[0] == '-' and a[1][0] == 'v' and a[2][0] == 'v':
        copyset.add(frozenset((a[1][1], a[2][1])))
ok = 0; tot = 0
for i, (Rw, c, Ex, aid, uc) in enumerate(conds):
    if info[i]['cls'] != 'diff' or info[i]['X'] != info[i]['Y']:
        continue
    tot += 1
    E = conds[i][2]
    if frozenset((E[1][1], E[2][1])) in copyset:
        ok += 1
print('(1) diff conditions with X==Y : %d ; backed by an explicit copy atom: %d' % (tot, ok))

# ---- (2) rank of the 3x2 (alpha,beta) matrix per block ----
byblk = defaultdict(list)
for (g, c, a1, a2, D1, D2) in rows:
    byblk[g].append((a1, a2))
bad = 0; mn = None; mx = 0
for g, ab in byblk.items():
    mins = [ab[i][0]*ab[j][1] - ab[j][0]*ab[i][1] for i, j in itertools.combinations(range(3), 2)]
    if any(m == 0 for m in mins):
        bad += 1
    m = max(abs(x) for x in mins); mx = max(mx, m)
    mn = m if mn is None else min(mn, m)
print('(2) blocks whose 3x2 (alpha,beta) matrix has a zero 2x2 minor: %d / %d' % (bad, len(byblk)))
print('    |minor| range: %d .. %d   (max < P: %s, < 2^48: %s)' % (mn, mx, mx < Pval, mx < 2**48))

# ---- (3) the gate family is a binary tree ----
fam = set()
for I, J in pure.values():
    fam.add(I); fam.add(J)
root = frozenset(range(256))
nodes = fam | {root}
parentok = 0; bij = Counter()
union2gate = {}
for g, (I, J) in pure.items():
    U = I | J
    if U in nodes:
        parentok += 1
    union2gate.setdefault(U, []).append(g)
print('(3) gates whose I∪J is itself a node of the family (or the root): %d / 255' % parentok)
print('    distinct unions: %d ; unions used by >1 gate: %d'
      % (len(union2gate), sum(1 for v in union2gate.values() if len(v) > 1)))
leaves = [x for x in nodes if len(x) == 1]
internal = [x for x in nodes if len(x) > 1]
print('    nodes: %d  (leaves %d, internal %d)  -- a binary tree on 256 leaves has 256+255=511'
      % (len(nodes), len(leaves), len(internal)))
print('    every internal node is the union of exactly one gate pair: %s'
      % (set(union2gate) == set(internal)))

# ---- (4) exact |S|-profile of each condition family ----
gate_of = {}
for i, (Rw, c, Ex, aid, uc) in enumerate(conds):
    if info[i]['cls'] in ('cong', 'offpin'):
        gate_of[i] = info[i]['gate']
merge = set(pure)
fams = defaultdict(lambda: defaultdict(int))
for i, (Rw, c, Ex, aid, uc) in enumerate(conds):
    d = info[i]; cls = d['cls']
    tag = '%s c%s' % (cls, '>1' if c > 1 else '=1')
    if cls in ('cong', 'offpin'):
        tag += ' [%s]' % ('merge' if d['gate'] in merge else 'dead')
    fams[tag]['n'] += 1
print('\n(4) condition census by family:')
for k in sorted(fams):
    print('    %-28s %d' % (k, fams[k]['n']))

# number of ACTIVE c>1 conditions as a function of S
cong_c = [info[i]['gate'] for i in range(len(conds))
          if info[i]['cls'] == 'cong' and conds[i][1] > 1 and info[i]['gate'] in merge]
offp_c = [info[i]['gate'] for i in range(len(conds))
          if info[i]['cls'] == 'offpin' and conds[i][1] > 1 and info[i]['gate'] in merge]
offp_dead = sum(1 for i in range(len(conds))
                if info[i]['cls'] == 'offpin' and conds[i][1] > 1 and info[i]['gate'] not in merge)
print('\n    c>1 congruences on merge blocks : %d' % len(cong_c))
print('    c>1 off-pins  on merge blocks    : %d' % len(offp_c))
print('    c>1 off-pins  on dead  blocks    : %d  (always active)' % offp_dead)

def liveset(S):
    return set(g for g, (I, J) in pure.items() if (I & S) and (J & S))
random.seed(7)
print('\n    |S|   live   act.cong(c>1)  act.offpin(c>1)  act.pin(c>1)  TOTAL active c>1')
for m in (1, 2, 4, 8, 16, 32, 64, 128, 192, 250, 256):
    for t in range(2):
        S = frozenset(random.sample(range(256), m))
        lv = liveset(S)
        ac = sum(1 for g in cong_c if g in lv)
        ao = sum(1 for g in offp_c if g not in lv) + offp_dead
        ap = m
        print('    %4d  %5d %12d %16d %13d %14d' % (m, len(lv), ac, ao, ap, ac + ao + ap + 11))
