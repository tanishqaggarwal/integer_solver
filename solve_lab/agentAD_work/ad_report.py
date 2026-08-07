"""Agent AD -- validate the extrapolation model on the measured points, then
extrapolate to n = 256 and say exactly how much the points can support.
"""
import json, os, sys
from math import comb, log, log2, exp


def refit_rho(clos, tot, min_tot=20, min_clos=5, min_pts=4):
    """Robust re-fit of the per-live-block factor from the stored counts.

    The raw fit in ad_scale.py uses every weight where the rate is in (0,1),
    which at n=8 can be five points with single-digit counts -- that is where the
    outlier rho values come from.  Here a weight only enters the fit if it has at
    least `min_tot` subsets and at least `min_clos` closing ones, and at least
    `min_pts` such weights must survive."""
    xs, ys = [], []
    for w in range(1, len(clos)):
        if tot[w] < min_tot or clos[w] < min_clos:
            continue
        r = clos[w] / tot[w]
        if r >= 1.0:
            continue
        xs.append(w - 1)
        ys.append(log(r))
    if len(xs) < min_pts:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return None
    return exp(sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den)

HERE = os.path.dirname(os.path.abspath(__file__))


def bmax_model(n, rho):
    best = -1
    for w in range(n + 1):
        if comb(n, w) * rho ** max(0, w - 1) >= 1.0:
            best = w
    return best


def rho_needed(n, B):
    """Largest rho with  C(n,w)*rho^(w-1) < 1  for every w > B."""
    lo, hi = 1e-9, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if bmax_model(n, mid) <= B:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    files = [f for f in os.listdir(HERE)
             if f.startswith('ad_scale') and f.endswith('.json')]
    rows = []
    for f in sorted(files):
        d = json.load(open(os.path.join(HERE, f)))
        for k, v in d.items():
            for rho, bm, cl, to in zip(v['rho'], v['bmax'], v['clos'], v['tot']):
                rows.append((k, v['n'], v['j0'], rho, bm, refit_rho(cl, to)))
    print('#### 1. MODEL VALIDATION  (model: #closing at weight w ~ C(n,w)*rho^(w-1))')
    print('     measured Bmax vs model Bmax, over every draw of every config')
    err = {}
    byn = {}
    for k, n, j0, rho, bm, rob in rows:
        if rho is None:
            continue
        pred = bmax_model(n, rho)
        e = pred - bm
        err[e] = err.get(e, 0) + 1
        byn.setdefault(n, []).append((rho, bm, pred, rob))
    tot = sum(err.values())
    for e in sorted(err):
        print('   model-Bmax - measured-Bmax = %+d : %5d draws (%.1f%%)'
              % (e, err[e], 100 * err[e] / tot))
    print('   |error| <= 1 in %.1f%% of %d draws'
          % (100 * sum(v for e, v in err.items() if abs(e) <= 1) / tot, tot))

    print('\n#### 2. rho BY SIZE  (is the per-live-block factor n-dependent?)')
    print('    n   draws   raw rho: min/med/max        ROBUST rho: min/med/max   n-Bmax')
    for n in sorted(byn):
        rs = sorted(r for r, _b, _p, _q in byn[n])
        rb = sorted(q for _r, _b, _p, q in byn[n] if q is not None)
        gaps = sorted(n - b for _r, b, _p, _q in byn[n])
        print('   %3d  %5d    %.4f %.4f %.4f   %s   %d..%d'
              % (n, len(rs), rs[0], rs[len(rs) // 2], rs[-1],
                 ('%.4f %.4f %.4f (%d)' % (rb[0], rb[len(rb) // 2], rb[-1], len(rb)))
                 if rb else 'n/a', gaps[0], gaps[-1]))

    print('\n#### 3. j=0 (CM by sqrt(-3), the real instance) vs generic curves')
    for j0 in (True, False):
        rs = sorted(r for _k, _n, j, r, _b, _q in rows if j == j0 and r is not None)
        if rs:
            print('   j0=%-5s draws=%4d  rho min=%.4f median=%.4f max=%.4f'
                  % (j0, len(rs), rs[0], rs[len(rs) // 2], rs[-1]))

    print('\n#### 4. EXTRAPOLATION TO n = 256')
    allr = sorted(r for _k, _n, _j, r, _b, _q in rows if r is not None)
    allrob = sorted(q for _k, _n, _j, _r, _b, q in rows if q is not None)
    for label, rho in (('most constraining RAW fit anywhere (n=8 small-sample tail)', allr[0]),
                       ('most constraining ROBUST fit anywhere', allrob[0]),
                       ('5th percentile, raw', allr[max(0, len(allr) // 20)]),
                       ('median, raw', allr[len(allr) // 2])):
        b = bmax_model(256, rho)
        print('   rho = %.4f  (%s)  ->  implied bound  w <= %d   [free bound is 255]'
              % (rho, label, b))
    print('\n   what rho would be needed at n = 256 for a USEFUL bound:')
    for B in (198, 148, 106, 56, 40, 24):
        r = rho_needed(256, B)
        print('     w <= %-4d needs rho <= %.4f   (per-live-block closure prob.)' % (B, r))
    print('\n   and the cost of that rho to the instance itself: the true subset,'
          '\n   of weight w, closes with probability rho^(w-1):')
    for B in (198, 148, 56, 24):
        r = rho_needed(256, B)
        print('     rho=%.4f -> a weight-%d true subset closes w.p. 2^%.1f'
              % (r, B, (B - 1) * log2(r)))


if __name__ == '__main__':
    main()
