"""lrun2 — fixes a REPORTING DEFECT in lrun.py and measures what the solves actually do.

`lattice.price` returned max(base_score, best_solution_score).  So its "39026 for all 127 row
subsets" said only *nothing beat the deliverable*; it hid whether the solve did anything at all
and what it cost.  A pricer that reports a floor instead of a measurement is the same failure
mode as reporting a maximum without its distribution, and I caught it in my own code.

This script reports, per row subset S:
  * whether the integer solve exists (it does for all 127 -- SNF says so),
  * whether re-propagating the solution ACTUALLY zeroes every row of S (verification, not trust),
  * the achieved score, unmaxed,
  * how many equations broke that were not in S.
"""
import sys, os, json, time, itertools, collections

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore, sparse                                    # noqa: E402
import lattice as LT                                           # noqa: E402

F = LT.FAILS_D
POOL, ATOMS = LT.reach_pool(F)
COLS = {}
DROP = []
for u in POOL:
    c = LT.column(u, LT.V_D, LT.BAD_D, F, LT.CM_D, LT.PIN)
    if c is None:
        DROP.append(u)
    elif c:
        COLS[u] = c
K = sorted(COLS)
RHS = {e: -LT.eqval(LT.BAD_D, e) for e in F}
print(f'frame D {LT.SCORE_D}/39033, rows {F}, knobs {len(K)} of pool {len(POOL)} '
      f'({len(DROP)} non-affine, dropped)', flush=True)


def solve_rows(S):
    rows = [{u: COLS[u][e] for u in K if e in COLS[u]} for e in S]
    rr = [RHS[e] for e in S]
    s, _, _ = sparse.solve_sparse(rows, rr, verbose=False, maxcore=400,
                                  maxcorebits=5_000_000)
    return s


out = []
t0 = time.time()
for r in range(1, len(F) + 1):
    for S in itertools.combinations(F, r):
        s = solve_rows(list(S))
        if s is None:
            out.append({'rows': list(S), 'solved': False})
            continue
        ch = {u: LT.V_D[u] + d for u, d in s.items() if d}
        bad, v = ieng.resid(LT.V_D, LT.BAD_D, ch, LT.PIN)
        sc = fscore.score(bad)
        fl = set(fscore.fails(bad))
        zeroed = [e for e in S if e not in fl]
        newbreak = sorted(fl - set(F))
        out.append({'rows': list(S), 'solved': True, 'score': sc,
                    'targets_zeroed': len(zeroed), 'targets': len(S),
                    'new_failures': len(newbreak), 'nknobs_moved': len(ch),
                    'still_failing': len(fl)})
        if sc > 39026:
            fn = f'M_lat2_{sc}_{"-".join(map(str, S))}.json'
            json.dump({f"x_{k}": int(v[k]) for k in range(LT.NV) if v[k] != 0}, open(fn, 'w'))
            print(f'  *** ABOVE 39026: {sc} buying {S} -> {fn} ***', flush=True)

print(f'\n{len(out)} row subsets, {time.time()-t0:.0f}s', flush=True)
okz = [o for o in out if o.get('solved') and o['targets_zeroed'] == o['targets']]
print(f'integer solve exists          : {sum(1 for o in out if o.get("solved"))} of {len(out)}',
      flush=True)
print(f'and re-propagation VERIFIES it : {len(okz)} of {len(out)}   '
      f'(every targeted row actually zero afterwards)', flush=True)

dist = collections.Counter(o['score'] for o in out if o.get('solved'))
print('\nachieved score (UNMAXED - this is what the solve really gets):', flush=True)
for s in sorted(dist, reverse=True):
    print(f'  {s}: {dist[s]}', flush=True)
print(f'\nbest achieved {max(dist)}  vs deliverable 39026', flush=True)

print('\nby number of rows bought:', flush=True)
for r in range(1, len(F) + 1):
    sub = [o for o in out if o.get('solved') and len(o['rows']) == r]
    print(f'  buy {r}: best {max(o["score"] for o in sub)}  '
          f'median new failures {sorted(o["new_failures"] for o in sub)[len(sub)//2]}',
          flush=True)

print('\nsingle-row prices (buy exactly one of the seven):', flush=True)
for o in out:
    if o.get('solved') and len(o['rows']) == 1:
        print(f'  eq {o["rows"][0]}: score {o["score"]}  failures {o["still_failing"]}  '
              f'new breaks {o["new_failures"]}  knobs moved {o["nknobs_moved"]}', flush=True)

json.dump({'knobs': K, 'dropped_nonaffine': DROP, 'rows': F, 'results': out},
          open('lrun2.json', 'w'), indent=1)
print('\nwrote lrun2.json', flush=True)
