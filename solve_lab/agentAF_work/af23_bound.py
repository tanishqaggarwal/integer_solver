#!/usr/bin/env python3
"""agent AF, step 23: (i) cross-check the deliverable's own live-count;
   (ii) the counterfactual ceiling: if a set Bad of merge blocks were unsatisfiable
        in live mode, what max |S| would follow?  (tree DP)"""
import sys, os, pickle, json
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import find
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
LAM = pickle.load(open(os.path.join(HERE, 'af_lam.pkl'), 'rb'))
T = pickle.load(open(os.path.join(HERE, 'af_tree.pkl'), 'rb'))
conds = C['conds']; info = M['info']; pure = LAM['pure']; sel = T['sel']
selidx = {find(s): i for i, s in enumerate(sel)}
merge = set(pure)

# ---- (i) the deliverable ----
d = json.load(open(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json')))
asg = {}
for k, v in d.items():
    asg[find(int(k[2:]) if k.startswith('x_') else int(k))] = int(v)
on = [i for s, i in selidx.items() if asg.get(s, 0) == 1]
tot = [i for s, i in selidx.items() if s in asg]
print('(i) deliverable: %d of 256 selectors assigned; ON = %s' % (len(tot), sorted(on)))
S = frozenset(on)
lv = [g for g, (I, J) in pure.items() if (I & S) and (J & S)]
print('    live merge blocks predicted by the gate law: %d   (|S|-1 = %d)' % (len(lv), max(len(S) - 1, 0)))

# ---- (ii) tree DP ----
# rebuild the tree: node = frozenset; children of the internal node I|J are I and J
children = {}
for g, (I, J) in pure.items():
    children[I | J] = (I, J, g)
root = frozenset(range(256))
assert root in children

def maxS(bad):
    memo = {}
    def f(node):
        if node in memo:
            return memo[node]
        if node not in children:
            memo[node] = 1; return 1
        I, J, g = children[node]
        a, b = f(I), f(J)
        memo[node] = max(a, b) if g in bad else a + b
        return memo[node]
    return f(root)

def cwith(pred):
    return set(info[i]['gate'] for i in range(len(conds))
               if conds[i][1] > 1 and pred(info[i]['cls'])
               and info[i].get('gate') in merge)

badcong = cwith(lambda c: c == 'cong')
badoff = cwith(lambda c: c == 'offpin')
print('\n(ii) counterfactual ceilings (tree DP over the recovered binary tree):')
print('     Bad = all 255 merge blocks            -> max |S| = %d' % maxS(merge))
print('     Bad = the 189 blocks with a c>1 congruence -> max |S| = %d' % maxS(badcong))
print('     Bad = the 121 blocks with a c>1 off-pin    -> max |S| = %d' % maxS(badoff))
print('     Bad = {root} only                     -> max |S| = %d' % maxS({children[root][2]}))
print('     Bad = {}                              -> max |S| = %d' % maxS(set()))
I, J, g = children[root]
print('     root split |I|=%d |J|=%d ; root gate carries a c>1 congruence: %s ; c>1 off-pin: %s'
      % (len(I), len(J), g in badcong, g in badoff))
known = [1, 2, 3, 5, 6, 7, 8, 17, 32, 64]
for nm, bad in (('189-cong', badcong), ('121-offpin', badoff), ('all-255', merge)):
    b = maxS(bad)
    bad_closures = [k for k in known if k > b]
    print('     ceiling %-11s = %3d   -> refuted by known closures %s'
          % (nm, b, bad_closures if bad_closures else 'NONE (consistent)'))
