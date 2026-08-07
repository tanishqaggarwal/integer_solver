"""Certificate: how few equations can ANY integer move of the full wide knob set leave failing?

Rows = the region R plus every outside equation any wide knob disturbs.  A move t leaves row i
satisfied iff (b + Mt)_i = 0, so the number of failing equations is the size of a minimum
"correction set": the complement of a maximum integrally solvable row subset.

Lower bound by disjoint minimal unsolvable subsets (MUSes).  Every unsolvable subset must lose at
least one of its rows, and disjoint MUSes must lose different rows, so

        min failing  >=  number of pairwise-disjoint MUSes.

Finding k >= 7 disjoint MUSes therefore PROVES no move of this knob set beats the deliverable's 7.
Each MUS is minimal by construction (deletion-based shrink, every step checked exactly), so the
certificate is re-checkable row by row.
"""
import json, time, sys
import zsolve
import optN
from optN import make, build, WIT
from widen import wide_knobs, build_wide


def shrink_to_mus(Z, S):
    """S unsolvable -> a minimal unsolvable subset of S (deletion-based)."""
    S = list(S)
    i = 0
    while i < len(S):
        T = S[:i] + S[i + 1:]
        if T and not Z.solvable(T):
            S = T
        else:
            i += 1
    return S


def find_mus(Z, pool):
    """Insertion-based: grow a subset until it is unsolvable, then shrink that (small) set.
    Far cheaper than deleting from the full row set, and the result is still a genuine MUS."""
    cur = []
    for r in pool:
        cur.append(r)
        if not Z.solvable(cur):
            return shrink_to_mus(Z, cur)
    return None


def disjoint_muses(Z, rows, limit=12, verbose=True):
    out = []
    pool = list(rows)
    while len(out) < limit:
        m = find_mus(Z, pool)
        if m is None:
            break
        out.append(sorted(m))
        if verbose:
            print('   MUS %d: size %d  rows %s' % (len(out), len(m), m), flush=True)
        pool = [r for r in pool if r not in set(m)]
        if not pool:
            break
    return out


def run(D, tag):
    st = make(list(D))
    d0 = build(st)
    Rl = d0['R']
    nR = len(Rl)
    knobs, outside = wide_knobs(st, Rl, verbose=False)
    rows, b, M, k = build_wide(st, Rl, knobs, outside)
    N = len(rows)
    print('\n=== %s ===  rows=%d (region %d + outside %d)  wide knobs=%d  base failing=%d'
          % (tag, N, nR, len(outside), k, len(st.fails)), flush=True)
    Z = zsolve.ZSolver(M, b, k)
    t0 = time.time()
    full = Z.solvable(range(N))
    print('  all %d rows simultaneously integrally solvable? %s' % (N, full), flush=True)
    if full:
        print('  *** the whole local model is closable -- construct and verify! ***', flush=True)
        return dict(D=list(D), rows=N, solvable=True)
    ms = disjoint_muses(Z, range(N))
    print('  disjoint MUSes found: %d  =>  min failing >= %d   (%.1fs)'
          % (len(ms), len(ms), time.time() - t0), flush=True)
    if len(ms) >= 7:
        print('  => CERTIFIED: no integer move of the %d wide knobs at this placement leaves '
              'fewer than 7 failing equations.' % k, flush=True)
    else:
        print('  => inconclusive from disjoint MUSes alone (bound %d < 7)' % len(ms), flush=True)
    return dict(D=list(D), rows=N, region=nR, outside=len(outside), knobs=k,
                nmus=len(ms), muses=ms, solvable=False,
                certified=len(ms) >= 7)


if __name__ == '__main__':
    res = []
    for D, tag in ((WIT, 'WITNESS %s' % WIT), ([28730], 'D=[28730]'), ([17499], 'D=[17499]'),
                   ([642, 28730, 31864], 'D=[642,28730,31864]'), ([], 'D=[] (all re-attached)')):
        try:
            res.append(run(D, tag))
        except Exception as e:
            print('  ERROR', e, flush=True)
    json.dump(res, open('runs/mus.json', 'w'), indent=1)
    print('\nwrote runs/mus.json')
