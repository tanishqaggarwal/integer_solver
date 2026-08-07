"""U20: price the cross-half ROUTE at every slot.

Construction (deliverable-style, honest pins):
  * turn on exactly two selectors  a in left subtree of beta,  b in right subtree
  * both leaf chains are made to CARRY leaf-`src`'s honest point  (the route lie)
  * both leaves keep their OWN honest pin constants (so no pin is violated, z=0)
  * beta's two inputs therefore coincide -> its residuals vanish -> its output is free
  * set beta's output so that the root carries the deliverable's own target
  * everything above beta is pass-through (sibling subtrees dead)
Score with checker.py's compiled equations through M's calibrated forward engine.
"""
import sys, pickle, time, json, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import umodel as U, ucrt as CR, uscore as SC, checker

DELIV = '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
v0 = checker.load_assignment(DELIV)
sd0 = SC.ENG.seed_of(v0)
# the 11 seed entries that are NOT produced by the build (target-drive machinery)
DRV = [642, 1329, 8731, 9118, 9413, 10903, 17325, 18956, 28730, 29854, 31864]
DRVSEED = {k: sd0[k] for k in DRV if k in sd0}
TGT = (v0[U.OUT[U.ROOT][0]['vab']], v0[U.OUT[U.ROOT][1]['vab']])


def build(beta, a, b, src):
    canon = U.LIFTC[src]
    rv = {a: canon, b: canon}
    bv = CR.betaval_for(beta, TGT)
    return U.assignment({a, b}, routeval=rv, beta=beta, betaval=bv)


def price(beta, a, b, src, drv=True):
    v, isl, valn = build(beta, a, b, src)
    s = SC.seed_of_build(v, DRVSEED if drv else None)
    n, vv = SC.score(s)
    return n, vv, s


def coincide(beta, a, b, src):
    """check the two projected inputs at beta really are equal (the whole point)"""
    canon = U.LIFTC[src]
    rv = {a: canon, b: canon}
    isl, valn = U.buildvals(set([a, b]), rv, None, None)
    la, lb = U.tree[beta]
    def proj(ch, side):
        pm = U.perm[(beta, side)]
        if valn[ch] is None or pm[0] is None or pm[1] is None:
            return None
        return (valn[ch][pm[0]], valn[ch][pm[1]])
    return proj(la, 'va'), proj(lb, 'vb')


if __name__ == '__main__':
    # ---- CONTROL: the deliverable's own slot/pair through the generalized function
    t0 = time.time()
    n, vv, s = price(U.ROOT, 24601, 2081, 24601)
    d = sum(1 for i in range(38748) if vv[i] != v0[i])
    print('CONTROL  beta=ROOT a=24601(e72) b=2081(e235) src=72 -> %d failing, %d vars differ (%.2fs)'
          % (n, d, time.time() - t0))
    PA, PB = coincide(U.ROOT, 24601, 2081, 24601)
    print('         inputs coincide at beta: %s' % (PA == PB))
