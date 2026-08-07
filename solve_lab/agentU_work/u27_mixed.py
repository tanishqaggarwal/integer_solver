"""U27: the mixed-axis route (each leaf honest on one coordinate, lying on the other).

u23 measured only two of the four axis-combinations.  If the price is additive over the two
coordinate routes then min over ALL carried values = min over the four combinations, and the
two mixed ones can be cheaper than either pure one.  Test additivity, then sweep the mixed
combinations at every slot for the best-known leaf pairs.
"""
import sys, time, pickle, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S
import umodel as U, ucrt as CR, uscore as SC


def pv(beta, a, b, canon):
    rv = {a: canon, b: canon}
    v, isl, valn = U.assignment({a, b}, routeval=rv, beta=beta,
                                betaval=CR.betaval_for(beta, S.TGT))
    return SC.score(SC.seed_of_build(v, S.DRVSEED))[0]


def four(beta, a, b):
    A, B = U.LIFTC[a], U.LIFTC[b]
    return (pv(beta, a, b, A),                              # b lies on both axes
            pv(beta, a, b, B),                              # a lies on both axes
            pv(beta, a, b, {'X': A['X'], 'Y': B['Y']}),     # a lies on Y, b lies on X
            pv(beta, a, b, {'X': B['X'], 'Y': A['Y']}))     # a lies on X, b lies on Y


if __name__ == '__main__':
    print('=== additivity test: does  pure_a + pure_b == mixed1 + mixed2 ? ===')
    tests = [(U.ROOT, 24601, 2081), (16102, 11368, 2081), (27994, 4287, 2081)]
    for beta in U.SLOTS[:6]:
        la, lb = U.tree[beta]
        tests.append((beta, sorted(U.LIVELEAF[la])[0], sorted(U.LIVELEAF[lb])[0]))
    for beta, a, b in tests:
        pa, pb, m1, m2 = four(beta, a, b)
        print('  beta=%-6d a=%-6d b=%-6d  pure_a=%-3d pure_b=%-3d mixed1=%-3d mixed2=%-3d'
              '  |  pa+pb=%-3d m1+m2=%-3d  additive=%s  min4=%d'
              % (beta, a, b, pa, pb, m1, m2, pa + pb, m1 + m2, pa + pb == m1 + m2,
                 min(pa, pb, m1, m2)))
