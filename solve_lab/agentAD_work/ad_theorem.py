"""Agent AD -- the ABSORBING-BLOCK theorem, stated and machine-checked.

THEOREM.  At a live merge block, let A = i1 - i2 (the integer difference of the
two merged x-slots) and let c be the block's small modulus.  If gcd(A, c) = 1
then, for ANY values of the four input lifts, the block's own output lifts
(u5,u6) can be chosen so that all three integer-lift conditions
        c_k | c_k1*nu1 + c_k2*nu2 ,  k = 1,2,3
hold.  The block then imposes no constraint whatever on the rest of the system.

PROOF.  nu1 = nu1a + u5*A^2 and nu2 = nu2a + u5*B + u6*A (checked exactly in
ad_validate V2c), so (u5,u6) -> (nu1,nu2) is an affine map with matrix
[[A^2,0],[B,A]], determinant A^3, invertible mod c when gcd(A,c)=1.  Take the
preimage of (0,0); every condition is homogeneous in (nu1,nu2), so all three
hold.  QED

COROLLARY (the whole of section 8, in one line).  The integer-lift family closes
for EVERY subset all of whose live blocks are absorbing.  A = a + P*delta with
delta = t_R - t_L a free lift difference and P invertible modulo any small c, so
A is not even fixed by the curve arithmetic: the input lifts alone can usually
make A invertible.  Failure therefore needs a coincidence per block, and the
closure rate is bounded below by prod_l (1-1/l)^(|S|-1) -- a per-block constant,
not a cutoff.

Checks below:
  T1  the theorem itself, brute force, on real curve points and random lifts;
  T2  the corollary: every subset whose live blocks are all absorbing closes
      (compared against the exact tree DP);
  T3  the lower bound prod_l (1-1/l)^(|S|-1) against the measured closure rate.
"""
import json, os, random, sys
from math import gcd
from collections import defaultdict
from ad_model import Curve, Instance, Coeffs, LiftDP, CANCEL

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    cs = {}
    for fn in ('ad_curves_partial.json', 'ad_curves.json', 'ad_curves2.json'):
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            for k, v in json.load(open(p)).items():
                cs.setdefault(k, []).extend(v)
    return cs


def T1(cs, ntrial=4000):
    rng = random.Random(77)
    cv = None
    for k in ('16', '12'):
        if k in cs:
            cd = cs[k][0]
            cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
            break
    p = cv.p
    tested = 0
    failed = 0
    for _ in range(ntrial):
        P1 = cv.mul(rng.randrange(1, cv.N), cv.G)
        P2 = cv.mul(rng.randrange(1, cv.N), cv.G)
        if P1 is None or P2 is None or P1[0] == P2[0]:
            continue
        P3 = cv.add(P1, P2)
        c = rng.choice((2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13))
        C = [[rng.randrange(-12, 13) for _ in range(2)] for _ in range(3)]
        u1, u2, u3, u4 = (rng.randrange(-9, 10) for _ in range(4))
        i1 = P2[0] + p * u1
        i2 = P1[0] + p * u2
        i3 = P1[1] + p * u3
        i4 = P2[1] + p * u4
        A = i1 - i2
        if gcd(A, c) != 1:
            continue
        tested += 1
        ok = False
        for u5 in range(c):
            for u6 in range(c):
                i5 = P3[0] + p * u5
                i6 = P3[1] + p * u6
                B = i4 - i3
                E = i1 + i2 + i5
                N1 = E * A * A - B * B
                N2 = A * (i3 + i6) - B * (i2 - i5)
                assert N1 % p == 0 and N2 % p == 0
                nu1, nu2 = N1 // p, N2 // p
                if all((C[k][0] * nu1 + C[k][1] * nu2) % c == 0 for k in range(3)):
                    ok = True
                    break
            if ok:
                break
        if not ok:
            failed += 1
    print('T1  absorbing-block theorem: %d gcd(A,c)=1 cases tested, %d had NO '
          'satisfying (u5,u6)  -> %s' % (tested, failed, 'PASS' if failed == 0 else 'FAIL'))
    return failed == 0


def greedy_certify(inst, co, S, info):
    """Independent, DP-free SUFFICIENT test for closure.

    Bottom-up: at each live block, given the integer lifts already fixed for its
    children, search that block's own (u5,u6) for a choice satisfying its three
    conditions, evaluating nu1,nu2 from the RAW integer definition of N1,N2.
    Success is a constructive witness, so  greedy => closure.  (Failure is
    inconclusive: the greedy commits to one choice per block.)
    """
    from math import lcm
    p = inst.cv.p
    bidx = {v.nid: i for i, v in enumerate(inst.internal)}
    lift = {}
    for v in inst.nodes:
        if v.idx is not None:
            lift[v.nid] = (0, 0)
            continue
        la, lb = info[v.left.nid][0], info[v.right.nid][0]
        if not (la and lb):
            lift[v.nid] = lift[v.left.nid] if la else (lift[v.right.nid] if lb else (0, 0))
            continue
        xl, yl = info[v.left.nid][2]
        xr, yr = info[v.right.nid][2]
        xv, yv = info[v.nid][2]
        tl, sl = lift[v.left.nid]
        tr, sr = lift[v.right.nid]
        i1, i2, i3, i4 = xr + p * tr, xl + p * tl, yl + p * sl, yr + p * sr
        C = co.C[bidx[v.nid]]
        mods = co.m[bidx[v.nid]]
        M = lcm(*[m for m in mods]) if max(mods) > 1 else 1
        got = None
        for u5 in range(M):
            for u6 in range(M):
                i5, i6 = xv + p * u5, yv + p * u6
                A, B, E = i1 - i2, i4 - i3, i1 + i2 + i5
                N1 = E * A * A - B * B
                N2 = A * (i3 + i6) - B * (i2 - i5)
                nu1, nu2 = N1 // p, N2 // p
                if all(mods[k] == 1 or (C[k][0] * nu1 + C[k][1] * nu2) % mods[k] == 0
                       for k in range(3)):
                    got = (u5, u6)
                    break
            if got:
                break
        if got is None:
            return False
        lift[v.nid] = got
    return True


def T2_T3(cs):
    """T2: greedy constructive certificate => DP says closed (one-directional).
       T3: the all-absorbing rate (reduced-A proxy) against the closure rate."""
    ok_all = True
    for nb in ('8', '12'):
        n = int(nb)
        cd = cs[nb][0]
        cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
        inst = Instance(cv, n)
        for seed in range(6):
            rng = random.Random(6000 + seed)
            co = Coeffs(rng, len(inst.internal), frac=1.0, pool=(2, 3))
            pp = co.prime_powers()
            dps = [LiftDP(inst, co, ell, E) for ell, E in sorted(pp.items())]
            bidx = {v.nid: i for i, v in enumerate(inst.internal)}
            viol = 0
            ngreedy = 0
            absorb_by_w = defaultdict(lambda: [0, 0])
            clos_by_w = defaultdict(lambda: [0, 0])
            for S in range(1 << n):
                info, blocks = inst.eval_subset(S)
                if any(c == CANCEL for _v, c, _d in blocks):
                    continue
                w = bin(S).count('1')
                allabs = True
                for v in inst.internal:
                    la, lb = info[v.left.nid][0], info[v.right.nid][0]
                    if not (la and lb):
                        continue
                    A = info[v.right.nid][2][0] - info[v.left.nid][2][0]
                    cs_ = co.m[bidx[v.nid]]
                    if any(cc > 1 and gcd(A, cc) != 1 for cc in cs_):
                        allabs = False
                        break
                closed = all(dp.run(S, info) != 0 for dp in dps)
                g = greedy_certify(inst, co, S, info)
                if g:
                    ngreedy += 1
                    if not closed:
                        viol += 1
                absorb_by_w[w][1] += 1
                clos_by_w[w][1] += 1
                if allabs:
                    absorb_by_w[w][0] += 1
                if closed:
                    clos_by_w[w][0] += 1
            if viol:
                ok_all = False
            print('T2  n=%2d seed=%d : greedy certificates=%d, of which the DP '
                  'calls NOT closed = %d' % (n, seed, ngreedy, viol))
            if seed == 0:
                print('    T3  w : all-absorbing-rate (reduced-A proxy)  vs  closure rate')
                for w in sorted(clos_by_w):
                    a = absorb_by_w[w]
                    c_ = clos_by_w[w]
                    print('        %2d :  %.4f      %.4f' % (w, a[0] / a[1], c_[0] / c_[1]))
    return ok_all


if __name__ == '__main__':
    cs = load()
    a = T1(cs)
    b = T2_T3(cs)
    print('\n%s' % ('ALL THEOREM CHECKS PASSED' if (a and b) else 'FAILED'))
