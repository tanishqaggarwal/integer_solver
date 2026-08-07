#!/usr/bin/env python3
"""agent AF, step 24: evaluate all 3707 lift conditions on the fleet's existing closures.
   This is a measurement on artefacts that already exist -- not a probe."""
import sys, os, pickle, json
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, Pval
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
LAM = pickle.load(open(os.path.join(HERE, 'af_lam.pkl'), 'rb'))
T = pickle.load(open(os.path.join(HERE, 'af_tree.pkl'), 'rb'))
conds = C['conds']; info = M['info']; pure = LAM['pure']; merge = set(pure)
sel = [find(s) for s in T['sel']]; selidx = {s: i for i, s in enumerate(sel)}
Rof = C['Rof']

def ev(n, A):
    t = n[0]
    if t == 'c':
        return n[1]
    if t == 'v':
        r = find(n[1])
        if r in val:
            return val[r]
        return A.get(r, 0)
    if t == 'neg':
        return -ev(n[1], A)
    a = ev(n[1], A); b = ev(n[2], A)
    return a + b if t == '+' else (a - b if t == '-' else a * b)

files = ['../best/new_instance_partial_39026.json'] + \
        ['../agentT_work/' + f for f in
         ['close_T2ctl.json', 'close_T3.json', 'close_T5.json', 'close_T6.json', 'close_T7.json',
          'close_T8.json', 'close_T17w.json', 'close_T32g.json', 'close_T64.json',
          'close_T128s59.json', 'close_T192s47.json', 'close_T250s31.json']]
print('%-26s %5s %6s %8s %8s %8s' % ('file', '|S|', 'live', 'lift-viol', 'of which', 'c>1'))
for f in files:
    p = os.path.join(HERE, f)
    if not os.path.exists(p):
        print('%-26s MISSING' % os.path.basename(f)); continue
    d = json.load(open(p))
    A = {}
    for k, v in d.items():
        A[find(int(k[2:]) if k.startswith('x_') else int(k))] = int(v)
    on = frozenset(i for s, i in selidx.items() if A.get(s, 0) == 1)
    lv = sum(1 for g, (I, J) in pure.items() if (I & on) and (J & on))
    bad = 0; badc = 0
    for i, (Rw, c, Ex, aid, uc) in enumerate(conds):
        v = ev(Ex, A)
        if v % (c * Pval) != 0:
            bad += 1
            if c > 1:
                badc += 1
    print('%-26s %5d %6d %8d %8s %8d' % (os.path.basename(f), len(on), lv, bad, '', badc))
