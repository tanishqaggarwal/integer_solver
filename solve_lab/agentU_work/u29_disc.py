"""U29: what IS the 5-equation discount at (honest=24601, lying=2081, beta=ROOT)?
And is it a property of the honest leaf, of the pair, or of the DRV seed I am holding fixed?
"""
import sys, time, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S
import umodel as U, uscore as SC, checker

CODES = SC.CODES


def fails(beta, a, b, src, drv=True):
    v, isl, valn = S.build(beta, a, b, src)
    sd = SC.seed_of_build(v, S.DRVSEED if drv else None)
    n, vv = SC.score(sd)
    return sorted(checker.evaluate_all(CODES, vv))


R = U.ROOT
f7 = fails(R, 24601, 2081, 24601)
f12 = fails(R, 47, 2081, 47)
print('honest=24601 (deliverable): %d failing  %s' % (len(f7), f7))
print('honest=47    (generic)    : %d failing  %s' % (len(f12), f12))
print('generic \\ deliverable = %s' % sorted(set(f12) - set(f7)))
print('deliverable \\ generic = %s' % sorted(set(f7) - set(f12)))

print('\n=== is the discount a property of the HONEST leaf 24601, or of the PAIR (72,235)? ===')
lb = U.tree[R][1]
B = sorted(U.LIVELEAF[lb])
prof = {}
t0 = time.time()
for b in B:
    prof[b] = len(fails(R, 24601, b, 24601))
c = collections.Counter(prof.values())
print('honest FIXED = 24601, lying varies over all %d leaves of the 78-half:' % len(B))
print('  distribution:', sorted(c.items()))
mn = min(prof.values())
print('  min %d at lying leaves %s' % (mn, [k for k, v in prof.items() if v == mn]))
print('  (%.0fs)' % (time.time() - t0))
pickle.dump(prof, open('u_disc_root.pkl', 'wb'))

print('\n=== does the discount survive dropping the deliverable DRV seed? ===')
for lab, drv in (('with DRV', True), ('no DRV', False)):
    a7 = len(fails(R, 24601, 2081, 24601, drv))
    a12 = len(fails(R, 47, 2081, 47, drv))
    print('  %-9s honest=24601 -> %-3d ; honest=47 -> %-3d ; discount %d'
          % (lab, a7, a12, a12 - a7))
