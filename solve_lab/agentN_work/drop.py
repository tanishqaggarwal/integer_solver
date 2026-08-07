"""How fast does the region open up as collateral is allowed?

Local model rows = region R (12) + outside equations any wide knob disturbs (139, all currently
zero).  Let W be the outside rows we allow to break.  Admissible moves are ker_Z(M_{outside\W});
within that lattice let g(W) = max number of integrally zeroable region rows.  Then

        failing = (|R| - g(W)) + |W|      and beating 7 needs   g(W) >= (|R|-6) + |W|.

Measured baseline g(empty) = 5 with |R| = 12, so |W| = 0 needs 6, |W| = 1 needs 7, and so on:
each row of collateral bought must buy strictly more than one region row.  This sweeps |W| = 0, 1,
2 exhaustively and reports the maximum g attained at each level.
"""
import json, time, itertools, sys
from flint import fmpz_mat
import zsolve
import optN
from optN import make, build, WIT
from widen import wide_knobs, build_wide, int_kernel

MAXD = int(sys.argv[1]) if len(sys.argv) > 1 else 2


def prep(D):
    st = make(list(D))
    d0 = build(st)
    Rl = d0['R']
    nR = len(Rl)
    knobs, outside = wide_knobs(st, Rl, verbose=False)
    rows, b, M, k = build_wide(st, Rl, knobs, outside)
    Mreg = [M[i] for i in range(nR)]
    breg = b[:nR]
    Mout = [M[i] for i in range(nR, len(rows))]
    return dict(nR=nR, k=k, Mreg=Mreg, breg=breg, Mout=Mout, nO=len(outside), st=st)


def g_of(P, keep):
    Mreg, breg, Mout, k, nR = P['Mreg'], P['breg'], P['Mout'], P['k'], P['nR']
    C = [Mout[i] for i in keep]
    K = int_kernel(C) if C else None
    if K is None:
        K = [[1 if a == bb else 0 for a in range(k)] for bb in range(k)]
    if not K:
        return 0, 0
    Proj = [[sum(Mreg[i][j] * v[j] for j in range(k)) for v in K] for i in range(nR)]
    opt, rws, exh, tests = zsolve.max_zero_rows(Proj, breg, len(K), nR)
    return opt, len(K)


def run(D, tag):
    P = prep(D)
    nO, nR = P['nO'], P['nR']
    need0 = nR - 6
    print('\n=== %s ===  |R|=%d outside=%d wide knobs=%d' % (tag, nR, nO, P['k']), flush=True)
    allout = list(range(nO))
    t0 = time.time()
    g0, dim0 = g_of(P, allout)
    print('  |W|=0: g=%d (kernel dim %d)  need %d  -> %s  (%.1fs)'
          % (g0, dim0, need0, 'BEATS' if g0 >= need0 else 'no', time.time() - t0), flush=True)
    best_overall = (nR - g0) + 0
    for d in range(1, MAXD + 1):
        need = need0 + d
        t0 = time.time()
        bg = -1
        bw = None
        n = 0
        for W in itertools.combinations(allout, d):
            keep = [i for i in allout if i not in set(W)]
            g, dim = g_of(P, keep)
            n += 1
            if g > bg:
                bg, bw = g, W
            fail = (nR - g) + d
            if fail < best_overall:
                best_overall = fail
            if g >= need:
                print('  *** |W|=%d W=%s gives g=%d >= %d : score %d ***'
                      % (d, list(W), g, need, 39033 - fail), flush=True)
        print('  |W|=%d: %d subsets, max g=%d (at W=%s), need %d -> %s ; best failing so far %d'
              ' (%.1fs)' % (d, n, bg, list(bw) if bw else None, need,
                            'BEATS' if bg >= need else 'no', best_overall, time.time() - t0),
              flush=True)
    print('  => over |W| <= %d the local model cannot do better than %d failing equations.'
          % (MAXD, best_overall), flush=True)
    return dict(D=list(D), nR=nR, nO=nO, g0=g0, maxd=MAXD, best_failing=best_overall)


if __name__ == '__main__':
    res = []
    for D, tag in ((WIT, 'WITNESS %s' % WIT), ([17499], 'D=[17499]')):
        res.append(run(D, tag))
    json.dump(res, open('runs/drop.json', 'w'), indent=1)
    print('\nwrote runs/drop.json')
