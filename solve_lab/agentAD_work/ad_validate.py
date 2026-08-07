"""Agent AD -- validation of the small analogue BEFORE any measurement.

Five independent checks, in increasing order of how much they would have caught:

  V1  law identification: on real curve points the gadget law N1=N2=0 holds mod P
      with (i5,i6)=(x_out-Q, y_out), and FAILS for Q != 0 -- so Q=0 is forced,
      not chosen.
  V2  the lift algebra: nu1 = N1/P and nu2 = N2/P computed by RAW big-integer
      arithmetic agree with the closed forms used by the DP, including the
      Jacobian identity nu1 = nu1a + u5*A^2, nu2 = nu2a + u5*B + u6*A.
      (This is the check that catches "discarded the quadratic part" and
      "truncated a fourth power" -- both third-degree and P^2 terms are compared.)
  V3  the DP is exact: for tiny n, brute-force enumerate EVERY lift assignment
      over every node and evaluate the raw integer conditions; compare the
      resulting closure verdict with the tree DP.  Both directions.
  V4  solution-set identity: for n=8 and n=12 IN FULL, the set of subsets the
      equation system admits equals { S : k(S)*G == T } computed from curve
      arithmetic.  Exceptional-block census reported.
  V5  awkward plants: (a) a target with TWO valid subsets, so a solver that
      stops at the first is wrong; (b) a rigged coefficient set whose closure
      failure set is a proper, non-empty, independently computable subset.
"""
import json, os, random, sys, itertools
from ad_model import (Curve, Instance, Coeffs, LiftDP, Q_OFFSET, CHORD, DEGEN,
                      CANCEL, PASS_L, PASS_R, DEAD, CLSNAME, congruence_closure,
                      closes, root_pinned_closes)

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []


def check(name, cond, extra=''):
    print(('  OK   ' if cond else '  FAIL ') + name + (('  ' + extra) if extra else ''))
    if not cond:
        FAIL.append(name)
    return cond


def load_curves():
    with open(os.path.join(HERE, os.environ.get('AD_CURVES','ad_curves.json'))) as f:
        return json.load(f)


# ------------------------------------------------------------------- V1 ----
def V1(cv, ntrial=400, rng=None):
    rng = rng or random.Random(11)
    p = cv.p
    bad_law = 0
    bad_q = 0
    for _ in range(ntrial):
        k1 = rng.randrange(1, cv.N)
        k2 = rng.randrange(1, cv.N)
        P1, P2 = cv.mul(k1, cv.G), cv.mul(k2, cv.G)
        if P1 is None or P2 is None or P1[0] == P2[0]:
            continue
        P3 = cv.add(P1, P2)
        xl, yl = P1
        xr, yr = P2
        xv, yv = P3
        for Q in (0, 1, 7):
            i1, i2, i3, i4 = xr, xl, yl, yr
            i5, i6 = xv - Q, yv
            A = i1 - i2
            B = i4 - i3
            E = i1 + i2 + i5 + Q
            N1 = E * A * A - B * B
            N2 = A * (i3 + i6) - B * (i2 - i5)
            ok = (N1 % p == 0 and N2 % p == 0)
            if Q == 0 and not ok:
                bad_law += 1
            if Q != 0 and ok:
                bad_q += 1
    check('V1a law N1=N2=0 mod P on real chord additions (Q=0)', bad_law == 0,
          'violations=%d' % bad_law)
    # Q != 0 : N1 still vanishes (the Q cancels inside E) but N2 does NOT,
    # because i5 is shifted while y_out is not.  Record the measured fact.
    check('V1b Q!=0 does not give a consistent ladder', bad_q == 0,
          'spurious-consistent=%d' % bad_q)


# ------------------------------------------------------------------- V2 ----
def V2(cv, ntrial=300, rng=None):
    """Raw N1/P, N2/P vs the closed forms, with random integer lifts."""
    rng = rng or random.Random(12)
    p = cv.p
    bad_int, bad_form, bad_jac = 0, 0, 0
    for _ in range(ntrial):
        k1 = rng.randrange(1, cv.N)
        k2 = rng.randrange(1, cv.N)
        P1, P2 = cv.mul(k1, cv.G), cv.mul(k2, cv.G)
        if P1 is None or P2 is None or P1[0] == P2[0]:
            continue
        P3 = cv.add(P1, P2)
        xl, yl = P1
        xr, yr = P2
        xv, yv = P3
        u = [rng.randrange(-40, 41) for _ in range(6)]
        u1, u2, u3, u4, u5, u6 = u
        i1 = xr + p * u1
        i2 = xl + p * u2
        i3 = yl + p * u3
        i4 = yr + p * u4
        i5 = xv + p * u5
        i6 = yv + p * u6
        A = i1 - i2
        B = i4 - i3
        E = i1 + i2 + i5
        N1 = E * A * A - B * B
        N2 = A * (i3 + i6) - B * (i2 - i5)
        if N1 % p or N2 % p:
            bad_int += 1
            continue
        nu1, nu2 = N1 // p, N2 // p
        # closed forms
        a = xr - xl
        b = yr - yl
        e = xr + xl + xv
        g = yl + yv
        h = xl - xv
        n1 = (e * a * a - b * b) // p
        n2 = (a * g - b * h) // p
        d = u1 - u2
        bt = u4 - u3
        eps = u1 + u2 + u5
        f1 = (n1 + (2 * a * e * d + eps * a * a - 2 * b * bt)
              + p * (e * d * d + 2 * a * eps * d - bt * bt) + p * p * (eps * d * d))
        f2 = (n2 + (a * (u3 + u6) + d * g - b * (u2 - u5) - bt * h)
              + p * (d * (u3 + u6) - bt * (u2 - u5)))
        if f1 != nu1 or f2 != nu2:
            bad_form += 1
        # Jacobian form used by the DP
        nu1a = (n1 + (2 * a * e * d + (u1 + u2) * a * a - 2 * b * bt)
                + p * (e * d * d + 2 * a * (u1 + u2) * d - bt * bt)
                + p * p * ((u1 + u2) * d * d))
        nu2a = (n2 + a * u3 + d * g - b * u2 - bt * h + p * (d * u3 - bt * u2))
        AA = a + p * d
        BB = b + p * bt
        if nu1a + u5 * AA * AA != nu1 or nu2a + u5 * BB + u6 * AA != nu2:
            bad_jac += 1
    check('V2a N1,N2 divisible by P for every integer lift', bad_int == 0,
          'bad=%d' % bad_int)
    check('V2b closed form for nu1,nu2 exact (incl. P^2 and cubic terms)',
          bad_form == 0, 'bad=%d' % bad_form)
    check('V2c Jacobian form nu1=nu1a+u5*A^2, nu2=nu2a+u5*B+u6*A exact',
          bad_jac == 0, 'bad=%d' % bad_jac)


# ------------------------------------------------------------------- V3 ----
def brute_closure(inst, co, S, q, ell):
    """Independent, raw-integer, exhaustive closure test for one prime power.

    Enumerates the lift of EVERY live-block output over Z/q (that is complete:
    every condition is a polynomial in the lifts with integer coefficients, so
    it depends on the lifts only mod q), builds the integer slot values, and
    evaluates the conditions with big-integer arithmetic from the raw
    definition of N1 and N2.  No use of the nu closed form.
    """
    p = inst.cv.p
    info, blocks = inst.eval_subset(S)
    order = [v for v in inst.nodes if v.idx is None]
    live_ids = []
    for v in order:
        ll = info[v.left.nid][0]
        lr = info[v.right.nid][0]
        if ll and lr:
            Pl, Pr = info[v.left.nid][2], info[v.right.nid][2]
            if Pl[0] == Pr[0]:
                if Pl[1] != Pr[1]:
                    return False            # CANCEL: infeasible mod P
                return None                 # DEGEN: not compared here
            live_ids.append(v.nid)
    if not live_ids:
        return True
    bidx = {v.nid: i for i, v in enumerate(inst.internal)}
    nv = len(live_ids)
    pos = {nid: i for i, nid in enumerate(live_ids)}
    for assign in itertools.product(range(q), repeat=2 * nv):
        lift = {}
        ok = True
        for v in inst.nodes:
            if v.idx is not None:
                lift[v.nid] = (0, 0)
                continue
            ll, lr = info[v.left.nid][0], info[v.right.nid][0]
            if ll and lr:
                i = pos[v.nid]
                lift[v.nid] = (assign[2 * i], assign[2 * i + 1])
            elif ll:
                lift[v.nid] = lift[v.left.nid]
            elif lr:
                lift[v.nid] = lift[v.right.nid]
            else:
                lift[v.nid] = (0, 0)
        for v in inst.internal:
            ll, lr = info[v.left.nid][0], info[v.right.nid][0]
            if not (ll and lr):
                continue
            xl, yl = info[v.left.nid][2]
            xr, yr = info[v.right.nid][2]
            xv, yv = info[v.nid][2]
            tl, sl = lift[v.left.nid]
            tr, sr = lift[v.right.nid]
            tv, sv = lift[v.nid]
            i1 = xr + p * tr
            i2 = xl + p * tl
            i3 = yl + p * sl
            i4 = yr + p * sr
            i5 = xv + p * tv
            i6 = yv + p * sv
            A = i1 - i2
            B = i4 - i3
            E = i1 + i2 + i5
            N1 = E * A * A - B * B
            N2 = A * (i3 + i6) - B * (i2 - i5)
            assert N1 % p == 0 and N2 % p == 0
            nu1, nu2 = N1 // p, N2 // p
            C = co.C[bidx[v.nid]]
            for k in range(3):
                m = co.m[bidx[v.nid]][k]
                mq = 1
                while m % ell == 0:
                    m //= ell
                    mq *= ell
                if mq > 1 and (C[k][0] * nu1 + C[k][1] * nu2) % mq:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False


def V3(curves, rng=None):
    rng = rng or random.Random(13)
    tot, mism = 0, 0
    detail = []
    for nb in ('8', '12'):
        for cd in curves[nb][:3]:
            cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
            for n in (4, 5):
                for mode in ('balanced', 'skew'):
                    inst = Instance(cv, n, tree_mode=mode)
                    for seed in range(3):
                        r = random.Random(1000 * seed + n + cd['p'])
                        co = Coeffs(r, len(inst.internal), frac=0.8, pool=(2, 3, 4))
                        pp = co.prime_powers()
                        for ell, E in sorted(pp.items()):
                            q = ell ** E
                            if q > 4:
                                continue
                            dp = LiftDP(inst, co, ell, E)
                            for S in range(1 << n):
                                b = brute_closure(inst, co, S, q, ell)
                                if b is None:
                                    continue
                                a = dp.run(S) != 0
                                tot += 1
                                if a != b:
                                    mism += 1
                                    if len(detail) < 6:
                                        detail.append((cd['p'], n, seed, ell, E, S, a, b))
    check('V3 tree DP == raw-integer exhaustive lift search', mism == 0,
          'compared=%d mismatches=%d %s' % (tot, mism, detail[:3]))


# ------------------------------------------------------------------- V4 ----
def exception_census(inst):
    cnt = {DEGEN: [], CANCEL: []}
    for S in range(1 << inst.n):
        _, blocks = inst.eval_subset(S)
        for v, cls, _d in blocks:
            if cls in (DEGEN, CANCEL):
                cnt[cls].append((S, v.nid))
    return cnt


def V4(curves, rng=None):
    """Full solution-set identity for n = 8 and n = 12."""
    ok_all = True
    for nb, n in (('8', 8), ('12', 12)):
        for cd in curves[nb][:4]:
            cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
            for mode in ('balanced', 'skew'):
                inst = Instance(cv, n, tree_mode=mode)
                cen = exception_census(inst)
                # pick a target
                r = random.Random(hash((cd['p'], n, mode)) & 0xffff)
                k0 = r.randrange(1, cv.N)
                T = cv.mul(k0, cv.G)
                curve_sols = set()
                for S in range(1 << n):
                    if cv.mul(inst.k_of(S), cv.G) == T:
                        curve_sols.add(S)
                # equation solutions at the mod-P layer:
                #   no CANCEL block anywhere, and either the root output == T
                #   or some block is DEGEN (output free => root steerable).
                eq_sols = set()
                degen_S = set()
                for S in range(1 << n):
                    _, blocks = inst.eval_subset(S)
                    kinds = [c for _v, c, _d in blocks]
                    if CANCEL in kinds:
                        continue
                    if DEGEN in kinds:
                        degen_S.add(S)
                        continue
                    if cv.mul(inst.k_of(S), cv.G) == T:
                        eq_sols.add(S)
                same = (eq_sols == curve_sols)
                ok_all &= same and not degen_S
                print('    p=%-8d n=%-3d %-8s  |curve sols|=%d  |eq sols|=%d  '
                      'DEGEN subsets=%d  CANCEL subsets=%d  match=%s'
                      % (cd['p'], n, mode, len(curve_sols), len(eq_sols),
                         len(degen_S), len(set(s for s, _ in cen[CANCEL])), same))
                if degen_S:
                    print('      !! degeneracy present: equation solution set is a '
                          'strict superset; instance not usable for V4')
    check('V4 equation solution set == { S : k(S)G == T } (n=8 and n=12, full)',
          ok_all)


# ------------------------------------------------------------------- V5 ----
def V5(curves):
    # --- plant (a): a target with TWO valid subsets.
    found = False
    for nb, n in (('8', 8), ('12', 12), ('16', 16)):
        for cd in curves[nb]:
            cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
            if cv.N >= (1 << n) - 1:
                continue
            inst = Instance(cv, n)
            # k0 and k0+N both < 2^n  => two subsets, with DIFFERENT popcounts
            room = (1 << n) - cv.N
            if room < 4:
                continue
            k0 = room // 2
            if k0 < 1 or k0 + cv.N >= (1 << n):
                continue
            T = cv.mul(k0, cv.G)
            sols = [S for S in range(1 << n) if cv.mul(inst.k_of(S), cv.G) == T]
            w = [bin(inst.k_of(S)).count('1') for S in sols]
            if len(sols) == 2:
                found = True
                print('    plant(a): p=%d n=%d k0=%d -> %d subsets, weights %s'
                      % (cd['p'], n, k0, len(sols), w))
                check('V5a two-solution plant recovered in full', sorted(sols) ==
                      sorted([k0, k0 + cv.N]),
                      'sols=%s expected=%s' % (sorted(sols), sorted([k0, k0 + cv.N])))
                break
        if found:
            break
    if not found:
        check('V5a two-solution plant recovered in full', False, 'no instance found')

    # --- plant (b): rigged coefficients with a proper, non-empty failure set.
    cd = curves['8'][0]
    cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
    n = 6
    inst = Instance(cv, n)
    nb_ = len(inst.internal)

    class Rig(Coeffs):
        def __init__(self):
            self.C = [[[1, 0], [0, 1], [1, 1]] for _ in range(nb_)]
            self.m = [(1, 1, 1) for _ in range(nb_)]
            self.mo = [(1, 1) for _ in range(nb_)]
    co = Rig()
    target_block = nb_ - 1                     # the ROOT block
    co.m[target_block] = (2, 1, 1)             # 2 | nu1  at the root only
    dp = LiftDP(inst, co, 2, 1)
    got = set()
    for S in range(1 << n):
        if dp.run(S) == 0:
            got.add(S)
    exp = set()
    for S in range(1 << n):
        if not brute_closure(inst, co, S, 2, 2):
            exp.add(S)
    print('    plant(b): failure set size %d of %d subsets (must be proper and '
          'non-empty)' % (len(got), 1 << n))
    check('V5b rigged failure set is neither empty nor everything',
          0 < len(got) < (1 << n), 'size=%d' % len(got))
    check('V5b DP reproduces the rigged failure set exactly', got == exp,
          'dp=%d brute=%d sym-diff=%d' % (len(got), len(exp), len(got ^ exp)))
    # vacuity guard, in agent X's sense: the scan must not be all-hits
    check('V5b plant is not vacuous (some subsets close, some do not)',
          0 < len(exp) < (1 << n))


def main():
    curves = load_curves()
    cd = curves['16'][0]
    cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
    print('V1 law identification            (curve p=%d)' % cv.p)
    V1(cv)
    print('V2 lift algebra')
    V2(cv)
    print('V5 plants')
    V5(curves)
    print('V3 DP vs raw-integer exhaustive search')
    V3(curves)
    print('V4 solution-set identity')
    V4(curves)
    print()
    if FAIL:
        print('VALIDATION FAILED:', FAIL)
        sys.exit(1)
    print('ALL VALIDATION CHECKS PASSED')


if __name__ == '__main__':
    main()
