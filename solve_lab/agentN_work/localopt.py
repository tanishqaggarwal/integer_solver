"""Exact settlement of the local model: can ANY integer move of the full wide knob set leave
fewer than 7 failing equations?

Rows = region R (currently 7 nonzero) + every outside equation any wide knob disturbs (139,
currently zero).  Any move t gives  failing = (region rows left nonzero) + (outside rows broken).
Writing A = set of region rows zeroed and W = set of outside rows broken,

        failing = (|R| - |A|) + |W|   <  7      requires   |W| <= |A| - (|R| - 6).

For every A the set U_A = { outside i : A u {i} not integrally solvable } must be broken, so
|W| >= |U_A|.  If |U_A| > |A| - (|R|-6) for every A, no move of the wide knob set beats 7.
The bound is a certificate, not a search: it is a necessary condition on A, checked exactly.
"""
import json, time, itertools, sys
import zsolve
import optN
from optN import make, build, WIT
from widen import wide_knobs, build_wide


def settle(D, tag, verbose=True):
    st = make(list(D))
    d0 = build(st)
    Rl = d0['R']
    nR = len(Rl)
    knobs, outside = wide_knobs(st, Rl, verbose=False)
    rows, b, M, k = build_wide(st, Rl, knobs, outside)
    nO = len(outside)
    base_fail = len(st.fails)
    need_slack = nR - 6            # |W| <= |A| - need_slack
    Z = zsolve.ZSolver(M, b, k)
    reg = list(range(nR))
    out = list(range(nR, nR + nO))
    print('\n=== %s ===' % tag, flush=True)
    print('  |R|=%d  outside touched=%d  wide knobs=%d  base failing=%d'
          % (nR, nO, k, base_fail), flush=True)
    # individually-compatible outside rows are the only ones that can ever be kept
    solo_bad = [i for i in out if not Z.solvable([i])]
    print('  outside rows not even individually keepable: %d (expected 0, they are 0 now)'
          % len(solo_bad), flush=True)
    t0 = time.time()
    survivors = []
    tested = 0
    minA = max(0, need_slack)          # |A| must be >= need_slack for |W|>=0 to be possible
    for size in range(nR, minA - 1, -1):
        budget = size - need_slack
        if budget < 0:
            continue
        for A in itertools.combinations(reg, size):
            if not Z.solvable(A):
                continue
            tested += 1
            bad = 0
            for i in out:
                if not Z.solvable(list(A) + [i]):
                    bad += 1
                    if bad > budget:
                        break
            if bad <= budget:
                survivors.append((sorted(A), bad, budget))
                print('  SURVIVOR A=%s |U_A|=%d budget=%d' % (sorted(A), bad, budget), flush=True)
    print('  region subsets that are themselves integrally solvable and pass the |U_A| test: %d'
          % len(survivors), flush=True)
    print('  (solvable region subsets examined: %d)  %.1fs' % (tested, time.time() - t0), flush=True)
    if not survivors:
        print('  => NO move of the full wide knob set leaves fewer than 7 failing equations '
              'at this placement.', flush=True)
    return survivors, dict(nR=nR, nO=nO, k=k, tested=tested, nsurv=len(survivors))


if __name__ == '__main__':
    res = []
    targets = [(WIT, 'WITNESS %s' % WIT), ([28730], 'D=[28730]'), ([17499], 'D=[17499]'),
               ([642, 28730, 31864], 'D=[642,28730,31864]')]
    for D, tag in targets:
        s, info = settle(D, tag)
        info['D'] = list(D)
        info['survivors'] = [[a, bb, bu] for a, bb, bu in s]
        res.append(info)
    json.dump(res, open('runs/localopt.json', 'w'), indent=1)
    print('\nwrote runs/localopt.json')
