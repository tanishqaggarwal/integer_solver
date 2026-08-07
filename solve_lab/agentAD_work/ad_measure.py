"""Agent AD -- the measurement: closure rate of the integer-lift families vs |S|.

Run only after ad_validate.py passes.

For each (curve, n, tree shape, coefficient draw) we compute, EXHAUSTIVELY over
all 2^n subsets (or by stratified sampling at n=20):

   * whether the CONGRUENCE lift family closes  (tree DP, exact)
   * whether the OFF-PIN lift family closes     (structural, see below)
   * the number of live merge blocks             (closed-form expectation |S|-1)

and report closure rate as a function of |S|, plus the single decisive
statistic:

   Bmax(draw) = max{ |S| : some subset of that weight closes }

If Bmax = n in every draw at every size, the integer-lift discipline places no
upper bound on |S| and section 8 of UPPER_BOUND_MAP is dead.  If Bmax < n, the
scaling of n - Bmax (or Bmax/n) with n is what extrapolates to n = 256.

Usage:  python3 ad_measure.py <mode>
  modes: main | perprime | stress | leaffree | n20 | quick
"""
import json, os, random, sys, time
from collections import defaultdict
from ad_model import Curve, Instance, Coeffs, LiftDP, CANCEL, DEGEN

HERE = os.path.dirname(os.path.abspath(__file__))
CURVES = os.environ.get('AD_CURVES', 'ad_curves.json')


def load():
    with open(os.path.join(HERE, CURVES)) as f:
        return json.load(f)


def popcount(x):
    return bin(x).count('1')


def modp_infeasible(inst, S):
    """True iff some live block has A=0,B!=0 (children opposite points)."""
    _, blocks = inst.eval_subset(S)
    return any(c == CANCEL for _v, c, _d in blocks)


def live_blocks(inst, S):
    _, blocks = inst.eval_subset(S)
    return sum(1 for _v, c, _d in blocks if c in (0, DEGEN, CANCEL))


def run_one(cv, n, mode, seed, frac, pool, leaf_free=False, subsets=None,
            check_liveblocks=False):
    """Returns dict weight -> [n_closed, n_total, n_rootpinned] plus diagnostics."""
    inst = Instance(cv, n, tree_mode=mode)
    rng = random.Random(seed)
    co = Coeffs(rng, len(inst.internal), frac=frac, pool=pool)
    leaf_mod = None
    if leaf_free:
        leaf_mod = [(rng.choice(pool) if rng.random() < frac else 1,
                     rng.choice(pool) if rng.random() < frac else 1)
                    for _ in range(n)]
    pp = co.prime_powers()
    dps = [LiftDP(inst, co, ell, E, leaf_free, leaf_mod) for ell, E in sorted(pp.items())]
    tab = defaultdict(lambda: [0, 0, 0])
    ncancel = 0
    lb_bad = 0
    it = range(1 << n) if subsets is None else subsets
    for S in it:
        info, blocks = inst.eval_subset(S)
        if any(c == CANCEL for _v, c, _d in blocks):
            ncancel += 1
            continue                       # mod-P infeasible, not a lift failure
        w = popcount(S)
        if check_liveblocks:
            nl = sum(1 for _v, c, _d in blocks if c in (0, DEGEN))
            if nl != max(0, w - 1):
                lb_bad += 1
        ok = True
        pinned = True
        for dp in dps:
            m = dp.run(S, info)
            if m == 0:
                ok = False
                pinned = False
                break
            if not (m & 1):
                pinned = False
        e = tab[w]
        e[1] += 1
        if ok:
            e[0] += 1
        if pinned:
            e[2] += 1
    return dict(tab), dict(nprimes=len(pp), primes=sorted(pp.items()),
                           ncancel=ncancel, lb_bad=lb_bad,
                           nmod=sum(1 for t in co.m for x in t if x > 1),
                           ncond=3 * len(inst.internal))


def summarize(name, per_draw, n):
    """per_draw: list of (tab, diag).  Print the full curve, not a statistic."""
    agg = defaultdict(lambda: [0, 0, 0])
    bmaxes = []
    for tab, _d in per_draw:
        bm = -1
        for w, (c, t, pin) in tab.items():
            a = agg[w]
            a[0] += c
            a[1] += t
            a[2] += pin
            if c > 0:
                bm = max(bm, w)
        bmaxes.append(bm)
    print('  %s   n=%d   draws=%d' % (name, n, len(per_draw)))
    print('    w :  closed/total    rate     rootpinned-rate')
    for w in sorted(agg):
        c, t, pin = agg[w]
        print('   %3d : %8d/%-8d %7.4f   %7.4f' % (w, c, t, c / t, pin / t))
    for i, (tab, _d) in enumerate(per_draw):
        print('    draw %d: %s   Bmax=%d' % (
            i, ' '.join('%d:%d/%d' % (w, tab[w][0], tab[w][1]) for w in sorted(tab)),
            bmaxes[i]))
    print('    Bmax per draw: %s' % bmaxes)
    print('    Bmax == n in %d/%d draws' % (sum(1 for b in bmaxes if b == n), len(bmaxes)))
    return agg, bmaxes


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    curves = load()
    out = {}

    def block(tag, sizes, ndraw, frac, pool, treemodes=('balanced', 'skew'),
              leaf_free=False, ncurve=5, sample=None):
        res = {}
        for nb in sizes:
            n = int(nb)
            for ci, cd in enumerate(curves[nb][:ncurve]):
                cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
                for tm in treemodes:
                    per = []
                    t0 = time.time()
                    for s in range(ndraw):
                        subs = None
                        if sample:
                            rr = random.Random(9000 + s)
                            subs = []
                            for w in range(n + 1):
                                for _ in range(sample):
                                    idx = rr.sample(range(n), w)
                                    subs.append(sum(1 << i for i in idx))
                            subs = list(set(subs))
                        tab, diag = run_one(cv, n, tm, 1000 + s, frac, pool,
                                            leaf_free, subs,
                                            check_liveblocks=(s == 0))
                        per.append((tab, diag))
                    key = '%s n=%d p=%d %s j0=%s' % (tag, n, cd['p'], tm, cd['j0'])
                    print('\n== %s   (%.1fs, primes=%s, nontrivial moduli %d/%d, '
                          'live-block-count mismatches %d)'
                          % (key, time.time() - t0, per[0][1]['primes'],
                             per[0][1]['nmod'], per[0][1]['ncond'], per[0][1]['lb_bad']))
                    agg, bmaxes = summarize(tag, per, n)
                    res[key] = dict(agg={str(k): v for k, v in agg.items()},
                                    bmax=bmaxes, n=n, j0=cd['j0'], tree=tm,
                                    primes=[list(x) for x in per[0][1]['primes']])
                    sys.stdout.flush()
        return res

    if mode == 'quick':
        out['quick'] = block('quick', ('8',), 4, 0.25, (2, 3, 4, 5), ncurve=2)
    elif mode == 'main':
        out['main'] = block('main', ('8', '12'), 12, 0.25,
                            (2, 3, 4, 5, 7, 8, 9, 11, 13))
    elif mode == 'main16':
        out['main16'] = block('main16', ('16',), 6, 0.25,
                              (2, 3, 4, 5, 7, 8, 9, 11, 13))
    elif mode == 'perprime':
        for q in (2, 3, 5, 7):
            out['q%d' % q] = block('q=%d' % q, ('8', '12'), 6, 0.25, (q,), ncurve=3)
    elif mode == 'stress':
        out['stress'] = block('stress', ('8', '12', '16'), 6, 1.0,
                              (2, 3, 4, 5, 7, 8, 9), ncurve=3)
    elif mode == 'stress2':
        out['stress2'] = block('stress2', ('8', '12'), 6, 1.0, (2,), ncurve=3)
        out['stress3'] = block('stress3', ('8', '12'), 6, 1.0, (2, 4, 8), ncurve=3)
    elif mode == 'leaffree':
        out['leaffree'] = block('leaffree', ('8', '12'), 6, 0.25,
                                (2, 3, 4, 5, 7, 8, 9), leaf_free=True, ncurve=3)
    elif mode == 'n20':
        out['n20'] = block('n20', ('20',), 4, 0.25, (2, 3, 4, 5, 7, 8, 9),
                           treemodes=('balanced',), ncurve=3, sample=400)
    with open(os.path.join(HERE, 'ad_measure_%s.json' % mode), 'w') as f:
        json.dump(out, f, indent=1)
    print('\nwritten ad_measure_%s.json' % mode)


if __name__ == '__main__':
    main()
