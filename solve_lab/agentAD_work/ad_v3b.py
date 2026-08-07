"""Agent AD -- V3b: the U-shape re-derived from RAW INTEGER arithmetic.

The headline measurement (closure rate is U-shaped in |S|, minimum near n/2,
back to 1 at |S| = n) is the whole result, so it is re-derived here without the
tree DP and without the nu closed form: for n = 7 every subset and every lift
assignment mod q is enumerated and the conditions are evaluated straight from
   N1 = E*A^2 - B^2 ,  N2 = A*(i3+i6) - B*(i2-i5)
on big integers.  Both the per-subset verdict AND the per-weight rate are
compared with the DP.
"""
import json, os, random, sys
from collections import defaultdict
from ad_model import Curve, Instance, Coeffs, LiftDP, CANCEL
from ad_validate import brute_closure

HERE = os.path.dirname(os.path.abspath(__file__))
cs = json.load(open(os.path.join(HERE, os.environ.get('AD_CURVES', 'ad_curves.json'))))

n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
mism = 0
tot = 0
for nb in ('8', '12'):
    for cd in cs[nb][:2]:
        cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
        inst = Instance(cv, n)
        for seed in range(3):
            r = random.Random(500 + seed + cd['p'])
            co = Coeffs(r, len(inst.internal), frac=0.6, pool=(2, 3))
            pp = co.prime_powers()
            for ell, E in sorted(pp.items()):
                q = ell ** E
                if q > 3:
                    continue
                dp = LiftDP(inst, co, ell, E)
                bt = defaultdict(lambda: [0, 0])
                dt = defaultdict(lambda: [0, 0])
                for S in range(1 << n):
                    b = brute_closure(inst, co, S, q, ell)
                    if b is None:
                        continue
                    a = dp.run(S) != 0
                    tot += 1
                    if a != b:
                        mism += 1
                    w = bin(S).count('1')
                    bt[w][1] += 1
                    dt[w][1] += 1
                    bt[w][0] += b
                    dt[w][0] += a
                print('p=%d seed=%d q=%d  brute: %s' %
                      (cd['p'], seed, q,
                       ' '.join('%d:%d/%d' % (w, bt[w][0], bt[w][1]) for w in sorted(bt))))
                assert bt == dt or dict(bt) == dict(dt), (dict(bt), dict(dt))
print('\ncompared %d (subset,config) pairs, %d mismatches' % (tot, mism))
print('PASS' if mism == 0 else 'FAIL')
