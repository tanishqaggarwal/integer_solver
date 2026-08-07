"""U28: RULE 9.  u24's reduction (price depends only on the lying leaf) is refuted by my own
control: u25 reports 12 at the ROOT with lying leaf 2081, but the deliverable IS the ROOT with
lying leaf 2081 and scores 7.  Measure the honest-leaf dependence properly.
"""
import sys, time, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S
import umodel as U

beta = U.ROOT
la, lb = U.tree[beta]
A = sorted(U.LIVELEAF[la]); B = sorted(U.LIVELEAF[lb])
print('ROOT=%d  |I|=%d |J|=%d ; 24601 in I: %s ; 2081 in J: %s'
      % (beta, len(A), len(B), 24601 in A, 2081 in B))

t0 = time.time()
prof = {}
for a in A:
    prof[a] = S.price(beta, a, 2081, a)[0]        # lying leaf fixed = 2081, honest varies
c = collections.Counter(prof.values())
print('\nlying leaf FIXED = 2081, honest leaf varies over all %d leaves of the 178-half:' % len(A))
print('  distribution:', sorted(c.items()))
print('  min %d at honest leaves %s' % (min(prof.values()),
      [k for k, v in prof.items() if v == min(prof.values())][:12]))
print('  price at honest=24601 (the deliverable): %d' % prof[24601])
print('  price at honest=A[0]=%d (what u25 used): %d' % (A[0], prof[A[0]]))
print('  (%.0fs)' % (time.time() - t0))
pickle.dump(prof, open('u_honest_root.pkl', 'wb'))
