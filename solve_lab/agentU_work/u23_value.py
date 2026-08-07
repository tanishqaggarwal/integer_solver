"""U23: does the route price depend on the VALUE carried, or only on the route geometry?
Caution (i) in RESUME_U §18/21 says the common point need not be a curve point.  Test it.
"""
import sys, random, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S
import umodel as U, ucrt as CR, uscore as SC, checker

TGT = S.TGT
random.seed(11)


def price_val(beta, a, b, canon):
    rv = {a: canon, b: canon}
    bv = CR.betaval_for(beta, TGT)
    v, isl, valn = U.assignment({a, b}, routeval=rv, beta=beta, betaval=bv)
    s = SC.seed_of_build(v, S.DRVSEED)
    n, vv = SC.score(s)
    return n


def trial(beta, a, b, label):
    print('== beta=%d a=%d b=%d  (%s)' % (beta, a, b, label))
    cases = [('honest a', U.LIFTC[a]), ('honest b', U.LIFTC[b])]
    others = [s for s in U.LIFTC if s not in (a, b)]
    for s in random.sample(others, 3):
        cases.append(('honest leaf %d' % s, U.LIFTC[s]))
    for i in range(4):
        cases.append(('random 296-bit off-curve #%d' % i,
                      {'X': random.getrandbits(296), 'Y': random.getrandbits(296)}))
    cases.append(('small (3,5)', {'X': 3, 'Y': 5}))
    cases.append(('zero', {'X': 0, 'Y': 0}))
    cases.append(('a.X with b.Y', {'X': U.LIFTC[a]['X'], 'Y': U.LIFTC[b]['Y']}))
    for lab, c in cases:
        t = time.time()
        n = price_val(beta, a, b, c)
        print('   %-28s -> %3d failing  (%.2fs)' % (lab, n, time.time() - t))


trial(U.ROOT, 24601, 2081, 'the deliverable slot+pair')
# a mid-depth slot from the u21 table
trial(16102, 11368, 2081, 'depth-3 slot')
trial(27994, 4287, 2081, 'depth-1 slot')
