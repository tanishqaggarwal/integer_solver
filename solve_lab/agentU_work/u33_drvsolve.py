"""U33: SOLVE for the DRV knobs at a slot instead of copying the deliverable's.

The discount lives in a 9-atom cluster {23616,23617,23618,36659..36664} whose free knobs are
the 11 DRV variables.  Four of the five discount equations are plain linear combinations of
atoms; each is therefore a low-degree polynomial in any single knob once the rest of the state
is fixed.  Fit it at k = 0,1,2, take the exact integer root, VERIFY the fit at k = 3, apply,
re-propagate, and score with checker.  Greedy to a fixpoint, from several starts.
"""
import sys, math, time, json, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work/mirror')
import u20_sweep as S
import umodel as U, uscore as SC, checker
import harness as H

ENG = SC.ENG
DISC = [2554, 6816, 8124, 8680, 9421]
DRV = [642, 1329, 8731, 9118, 9413, 10903, 17325, 18956, 28730, 29854, 31864]
NS = {'__builtins__': {}}


def eqsum(v, e):
    ns = {'v': v, '__builtins__': {}}
    issq, outer, terms = H.eqt[e]
    s = 0
    for c, a in terms:
        s += c if a < 0 else c * eval(H.acodes[a], ns)
    return s


def score(seed):
    v = ENG.forward(seed)
    return len(checker.evaluate_all(SC.CODES, v)), v


def solve_knob(seed, e, d):
    """exact integer roots of  eqsum(e) == 0  as a function of knob d, with a degree check"""
    y = []
    for t in (0, 1, 2, 3):
        s2 = dict(seed); s2[d] = t
        y.append(eqsum(ENG.forward(s2), e))
    c0, c1, c2, c3 = y
    A2 = c2 - 2 * c1 + c0
    if A2 == 0:
        sl = c1 - c0
        if c0 + 3 * sl != c3:
            return []                      # not linear in d -> refuse
        if sl == 0 or c0 % sl:
            return []
        return [-c0 // sl]
    A = A2 // 2; B = c1 - c0 - A; C = c0
    if A * 9 + B * 3 + C != c3:
        return []                          # not quadratic in d -> refuse
    disc = B * B - 4 * A * C
    if disc < 0:
        return []
    r = math.isqrt(disc)
    if r * r != disc:
        return []
    return sorted({(-B + s) // (2 * A) for s in (r, -r) if (-B + s) % (2 * A) == 0})


def greedy(seed0, tag, rounds=6, verbose=True):
    seed = dict(seed0)
    best, _ = score(seed)
    if verbose:
        print('   [%s] start %d failing' % (tag, best))
    for it in range(rounds):
        moved = False
        for d in DRV:
            for e in DISC:
                for root in solve_knob(seed, e, d):
                    s2 = dict(seed)
                    if root == 0:
                        s2.pop(d, None)
                    else:
                        s2[d] = root
                    n, _ = score(s2)
                    if n < best:
                        best, seed, moved = n, s2, True
                        if verbose:
                            print('   [%s] knob x_%d := root of eq%d -> %d failing'
                                  % (tag, d, e, n))
        if not moved:
            break
    return best, seed


if __name__ == '__main__':
    v0 = checker.load_assignment(
        '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
    print('=== sanity: eqsum on the deliverable ===')
    for e in DISC:
        print('  eq%-6d sum=%s' % (e, eqsum(v0, e)))
    print('  (all zero => the deliverable satisfies all five, as checker says)')

    print('\n=== CONTROL: can the solver recover the discount at the ROOT from a GENERIC pair? ===')
    beta, a, b = U.ROOT, 47, 2081
    v, isl, valn = S.build(beta, a, b, a)
    for tag, extra in (('deliverable DRV', S.DRVSEED), ('no DRV', None)):
        sd = SC.seed_of_build(v, extra)
        n, sol = greedy(sd, tag)
        print('   [%s] FINAL %d failing' % (tag, n))
        if n <= 11:
            vv = ENG.forward(sol)
            json.dump({("x_%d" % i): vv[i] for i in range(38748) if vv[i] != 0},
                      open('u_drv_root47_%d.json' % n, 'w'))
            print('   dumped u_drv_root47_%d.json' % n)
