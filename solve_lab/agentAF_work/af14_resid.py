#!/usr/bin/env python3
"""agent AF, step 14: what is the residual wire in  c*P | L*B ?  Is it free, or computed?"""
import sys, os, pickle, json
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
conds = C['conds']; info = M['info']

def raw(n):
    if n[0] == 'c':
        return str(n[1]) if abs(n[1]) < 10**12 else 'BIG'
    if n[0] == 'v':
        return 'x%d' % find(n[1])
    if n[0] == 'neg':
        return '-' + raw(n[1])
    return '(%s %s %s)' % (raw(n[1]), n[0], raw(n[2]))

st = Counter()
freeB = 0
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    d = info[i]
    if d['cls'] not in ('cong', 'offpin'):
        continue
    B = d['other']
    dl = defc.get(B, [])
    st[(d['cls'], 'c>1' if c > 1 else 'c=1', len(dl))] += 1
    if not dl:
        freeB += 1
print('(class, c, #definitions of the residual wire) census:')
for k in sorted(st):
    print('   %-8s %-4s ndefs=%d : %d' % (k[0], k[1], k[2], st[k]))

print('\nresidual wire shapes (1 definition):')
sh = Counter()
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    d = info[i]
    if d['cls'] not in ('cong', 'offpin'):
        continue
    dl = defc.get(d['other'], [])
    if len(dl) == 1:
        sh[(d['cls'], shape_of(expand(('v', d['other']), 1)))] += 1
for k, v in sh.most_common(14):
    print('   %5d  %-8s %s' % (v, k[0], k[1]))

print('\nexamples, congruence with c>1 (residual expanded 3 levels):')
n = 0
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    d = info[i]
    if d['cls'] != 'cong' or c == 1:
        continue
    print('  c=%d   L=x%d   B=x%d  =  %s' % (c, d['gate'], d['other'], raw(expand(('v', d['other']), 3))[:230]))
    n += 1
    if n >= 6:
        break
print('\nexamples, off-pin with c>1:')
n = 0
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    d = info[i]
    if d['cls'] != 'offpin' or c == 1:
        continue
    print('  c=%d   1-L=x%d  B=x%d  =  %s' % (c, d['gate'], d['other'], raw(expand(('v', d['other']), 3))[:230]))
    n += 1
    if n >= 6:
        break
