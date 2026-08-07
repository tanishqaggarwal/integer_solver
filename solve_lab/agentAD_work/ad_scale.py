"""Agent AD -- how does the closure cutoff scale with n?

The measurement in ad_measure.py shows closure rate falling with |S|.  The only
thing that matters for section 8 of UPPER_BOUND_MAP is whether that fall is fast
enough, at n = 256, to empty out every weight above some useful B.  This script
measures the fall exactly, at n = 8,10,12,14,16,18,20, and tests an explicit
extrapolation model against those points before using it.

MODEL (stated before fitting, tested after):

    #closing subsets of weight w  ~  C(n,w) * rho^(w-1)

    i.e. each of the w-1 live merge blocks independently fails to absorb the
    integer-lift conditions with probability 1-rho.  (Live-block count = |S|-1
    exactly, checked in ad_measure.run_one.)

    => Bmax(n) = max{ w : C(n,w)*rho^(w-1) >= 1 }.

If the model reproduces the measured Bmax at every n, extrapolating it to 256 is
at least a defensible statement of what these points can support -- and the
extrapolation is then reported with the same caveat AB attached to d_reg.
"""
import json, os, random, sys, time
from collections import defaultdict
from math import comb, log
from ad_model import Curve, Instance, Coeffs, LiftDP, CANCEL

HERE = os.path.dirname(os.path.abspath(__file__))


def load_all():
    cs = {}
    for fn in ('ad_curves_partial.json', 'ad_curves.json', 'ad_curves2.json'):
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            for k, v in json.load(open(p)).items():
                cs.setdefault(k, [])
                have = {(c['p'], c['a'], c['b']) for c in cs[k]}
                for c in v:
                    if (c['p'], c['a'], c['b']) not in have:
                        cs[k].append(c)
    return cs


def counts_by_weight(cv, n, mode, seed, frac, pool):
    inst = Instance(cv, n, tree_mode=mode)
    rng = random.Random(seed)
    co = Coeffs(rng, len(inst.internal), frac=frac, pool=pool)
    pp = co.prime_powers()
    dps = [LiftDP(inst, co, ell, E) for ell, E in sorted(pp.items())]
    clos = [0] * (n + 1)
    tot = [0] * (n + 1)
    for S in range(1 << n):
        info, blocks = inst.eval_subset(S)
        if any(c == CANCEL for _v, c, _d in blocks):
            continue
        w = bin(S).count('1')
        tot[w] += 1
        ok = True
        for dp in dps:
            if dp.run(S, info) == 0:
                ok = False
                break
        if ok:
            clos[w] += 1
    return clos, tot, sorted(pp.items())


def fit_rho(clos, tot):
    """Least-squares slope of log(rate) against w over the range where the rate
    is measurable (>= 1/total and < 1)."""
    xs, ys = [], []
    for w in range(1, len(clos)):
        if tot[w] == 0 or clos[w] == 0:
            continue
        r = clos[w] / tot[w]
        if r >= 1.0:
            continue
        xs.append(w - 1)
        ys.append(log(r))
    if len(xs) < 3:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return None
    import math
    return math.exp(num / den)


def bmax_model(n, rho):
    best = -1
    for w in range(n + 1):
        if comb(n, w) * rho ** max(0, w - 1) >= 1.0:
            best = w
    return best


def main():
    cs = load_all()
    configs = []
    for spec in sys.argv[1:] or ['2:1.0', '2,3:1.0', '2,3,4,5,7,8,9:1.0',
                                 '2,3,4,5,7,8,9,11,13:0.25']:
        pool, frac = spec.split(':')
        configs.append((tuple(int(x) for x in pool.split(',')), float(frac)))
    sizes = [int(k) for k in sorted(cs, key=int) if int(k) <= 20]
    out = {}
    for pool, frac in configs:
        print('\n############ pool=%s frac=%s' % (pool, frac))
        for n in sizes:
            pick = []
            for want in (True, False):
                for c in cs[str(n)]:
                    if bool(c['j0']) == want:
                        pick.append(c)
                        break
            for cd in pick:
                cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
                ndraw = 8 if n <= 14 else (5 if n <= 18 else 3)
                rows = []
                t0 = time.time()
                for s in range(ndraw):
                    clos, tot, pp = counts_by_weight(cv, n, 'balanced',
                                                     4000 + s, frac, pool)
                    rho = fit_rho(clos, tot)
                    bm = max([w for w in range(n + 1) if clos[w] > 0])
                    rows.append((clos, tot, rho, bm, pp))
                print('\n-- n=%2d p=%d j0=%s   (%.1fs, primes=%s)'
                      % (n, cd['p'], cd['j0'], time.time() - t0, rows[0][4]))
                for i, (clos, tot, rho, bm, _pp) in enumerate(rows):
                    pred = bmax_model(n, rho) if rho else n
                    print('   draw %d  Bmax=%-3d rho=%-7s model-Bmax=%-3d  counts=%s'
                          % (i, bm, ('%.4f' % rho) if rho else 'n/a', pred,
                             ','.join(str(c) for c in clos)))
                rr = [r[2] for r in rows if r[2]]
                bms = [r[3] for r in rows]
                print('   rho: min=%s median=%s ; Bmax: %s ; n-Bmax: %s'
                      % (('%.4f' % min(rr)) if rr else 'n/a',
                         ('%.4f' % sorted(rr)[len(rr) // 2]) if rr else 'n/a',
                         bms, [n - b for b in bms]))
                out['pool%s_frac%s_n%d_p%d' % ('-'.join(map(str, pool)), frac, n, cd['p'])] = \
                    dict(n=n, p=cd['p'], j0=cd['j0'],
                         clos=[r[0] for r in rows], tot=[r[1] for r in rows],
                         rho=[r[2] for r in rows], bmax=bms)
                sys.stdout.flush()
    with open(os.path.join(HERE, os.environ.get('AD_OUT','ad_scale.json')), 'w') as f:
        json.dump(out, f, indent=1)
    print('\nwritten ad_scale.json')


if __name__ == '__main__':
    main()
