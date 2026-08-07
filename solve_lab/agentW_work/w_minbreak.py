"""TASK 1b: EXACT minimum-break per bought set, by implicit hitting set over unsat cores.

For a target P (subset of the 7 failing rows we insist on buying), define
    minbreak(P) = min |B|, B subset of SAT, such that (SAT \ B) + P has an INTEGER solution.
Gain = |P| - minbreak(P).  An improvement over 39,026 exists iff some P has gain > 0.

Brute force costs C(168,b) solves per P; this costs a few hundred and returns the EXACT
minimum instead of a yes/no at a cutoff.  Oracle is the same exact FLINT-HNF integer solve O
used (audited in w_oracle.py: every None is 'core infeasible', never a give-up), so the oracle
is MONOTONE: adding rows never turns infeasible into feasible.

Soundness: every infeasible kept-set C yields a core -- at least one row of C must be broken.
A B that hits no accumulated core cannot be a solution, so the minimum hitting set size is a
valid LOWER bound on minbreak; the loop stops when that bound is realised (exact answer) or
exceeds the cutoff (proved minbreak > cutoff).
"""
import sys, os, itertools, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S

INF = 10 ** 9
SATSET = set(S.SAT)


def feas_keep(P, keep, cache):
    """is (keep + P) integer-feasible?  keep is a frozenset of SAT row names."""
    key = (P, keep)
    r = cache.get(key)
    if r is None:
        r = S.solve(sorted(keep, key=str) + list(P)) is not None
        cache[key] = r
    return r


def min_hit(cores, cutoff):
    """smallest hitting set of `cores` with size <= cutoff, or None."""
    if not cores:
        return set()
    best = [None]
    for k in range(0, cutoff + 1):
        best[0] = None
        def rec(chosen):
            if best[0] is not None:
                return
            rem = [c for c in cores if not (c & chosen)]
            if not rem:
                best[0] = set(chosen); return
            if len(chosen) >= k:
                return
            for r in sorted(min(rem, key=len), key=str):
                chosen.add(r); rec(chosen); chosen.discard(r)
                if best[0] is not None: return
        rec(set())
        if best[0] is not None:
            return best[0]
    return None


def shrink(keep, P, cache):
    """minimal core: subset C of `keep` with C+P still infeasible."""
    C = set(keep)
    for r in sorted(keep, key=str):
        if not feas_keep(P, frozenset(C - {r}), cache):
            C.discard(r)
    return frozenset(C)


def minbreak(P, cutoff, log=None):
    cache = {}
    if S.solve(list(P)) is None:
        return (INF, None, {'note': 'P infeasible on its own'})
    cores, it = [], 0
    while True:
        B = min_hit(cores, cutoff)
        if B is None:
            return (None, None, {'cores': len(cores), 'coresizes': sorted(len(c) for c in cores),
                                 'solves': len(cache)})
        keep = frozenset(SATSET - B)
        if feas_keep(P, keep, cache):
            return (len(B), sorted(B, key=str), {'cores': len(cores), 'solves': len(cache)})
        c = shrink(keep, P, cache)
        cores.append(set(c)); it += 1
        if log:
            print('    core %d size %d = %s   (|B| tried %d, %d solves)'
                  % (it, len(c), sorted(c, key=str)[:12], len(B), len(cache)), file=log, flush=True)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    LOG = open(os.path.join(S.HERE, 'w_minbreak_%s.log' % which), 'w', buffering=1)
    def say(*a): print(*a, file=LOG, flush=True)
    sizes = [int(x) for x in which.split(',')] if which != 'all' else [1, 2, 3, 4, 5, 6, 7]
    res, t0 = {}, time.time()
    for k in sizes:
        for P in itertools.combinations(S.FAIL, k):
            tp = time.time()
            cutoff = k - 1
            v, B, st = minbreak(P, cutoff, log=LOG)
            if v is None:
                say('buy %-34s : minbreak > %d  -> NO GAIN  (%s, %.0fs)' % (list(P), cutoff, st, time.time()-tp))
                res[str(list(P))] = {'minbreak_gt': cutoff, **st}
            elif v == INF:
                say('buy %-34s : INFEASIBLE ALONE -> NO GAIN (%.0fs)' % (list(P), time.time()-tp))
                res[str(list(P))] = {'infeasible_alone': True}
            else:
                say('buy %-34s : minbreak = %d B=%s -> GAIN %d (%s, %.0fs)' % (list(P), v, B, k-v, st, time.time()-tp))
                res[str(list(P))] = {'minbreak': v, 'B': B, 'gain': k-v, **st}
                sol = S.solve(sorted(SATSET - set(B), key=str) + list(P))
                S.price(sol, 'realising buy %s break %s' % (list(P), B), log=LOG, tagfile='mb')
            json.dump(res, open(os.path.join(S.HERE, 'w_minbreak_%s.json' % which), 'w'), indent=1)
    say('TOTAL %.0fs' % (time.time()-t0)); say('DONE')
