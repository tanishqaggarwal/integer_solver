#!/usr/bin/env python3
"""agent AF, step 13: map every one of the 3707 conditions to its block and gate polarity."""
import sys, os, pickle
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import find, pp, atoms, defc
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
B = pickle.load(open(os.path.join(HERE, 'af_blocks.pkl'), 'rb'))
G = pickle.load(open(os.path.join(HERE, 'af_gates.pkl'), 'rb'))
LAM = pickle.load(open(os.path.join(HERE, 'af_lam.pkl'), 'rb'))
T = pickle.load(open(os.path.join(HERE, 'af_tree.pkl'), 'rb'))
conds = C['conds']; gates = B['gates']; prod = B['prod']; pins = B['pins']; diff = B['diff']
gf = G['gatefn']; pure = LAM['pure']; sel = T['sel']; selidx = {s: i for i, s in enumerate(sel)}
notg = {v: k for k, v in gates.items()}     # (1-L) var -> L var

# per-condition classification
info = {}                     # cond index -> dict
for (i, k, A, Bv) in prod:
    if A in gates or A in notg:
        g, other = A, Bv
    elif Bv in gates or Bv in notg:
        g, other = Bv, A
    else:
        info[i] = dict(cls='prod-nogate'); continue
    L = g if g in gates else notg[g]
    pol = 1 if g in gates else 0          # 1: gated by L (congruence) ; 0: gated by (1-L) (off-pin)
    info[i] = dict(cls='cong' if pol else 'offpin', gate=L, other=other)
for (i, k, s, w, Cst) in pins:
    info[i] = dict(cls='pin', selector=find(s))
for (i, k, X, Y, sg) in diff:
    info[i] = dict(cls='diff', X=X, Y=Y)
for i in range(len(conds)):
    info.setdefault(i, dict(cls='other'))
print('condition classes:', dict(Counter(v['cls'] for v in info.values())))

cbyc = defaultdict(Counter)
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    cbyc[info[i]['cls']][('c>1' if c > 1 else 'c=1')] += 1
for k in sorted(cbyc):
    print('  %-12s %s' % (k, dict(cbyc[k])))

# per-block breakdown
blk = defaultdict(lambda: dict(cong=[], offpin=[]))
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    d = info[i]
    if d['cls'] in ('cong', 'offpin'):
        blk[d['gate']][d['cls']].append(c)
print('\nblocks seen: %d' % len(blk))
print('per block (#cong, #offpin):', dict(Counter((len(v['cong']), len(v['offpin'])) for v in blk.values())))
print('per block #cong with c>1 :', dict(sorted(Counter(sum(1 for c in v['cong'] if c > 1) for v in blk.values()).items())))
print('per block #offpin with c>1:', dict(sorted(Counter(sum(1 for c in v['offpin'] if c > 1) for v in blk.values()).items())))

live = set(pure)                     # the 255 merge gates
dead = set(gates) - live
print('\nmerge (AND) gates: %d   permanently-dead gates: %d' % (len(live), len(dead)))
for name, gs in (('MERGE', live), ('CONST0', dead)):
    nc = sum(sum(1 for c in blk[g]['cong'] if c > 1) for g in gs)
    no = sum(sum(1 for c in blk[g]['offpin'] if c > 1) for g in gs)
    nc1 = sum(sum(1 for c in blk[g]['cong'] if c == 1) for g in gs)
    no1 = sum(sum(1 for c in blk[g]['offpin'] if c == 1) for g in gs)
    print('  %-7s  congruence c>1 %4d / c=1 %4d      off-pin c>1 %4d / c=1 %4d'
          % (name, nc, nc1, no, no1))

pickle.dump({'info': info, 'blk': dict(blk), 'live': live, 'dead': dead},
            open(os.path.join(HERE, 'af_map.pkl'), 'wb'), 2)
