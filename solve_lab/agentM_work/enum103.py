"""Enumerate 4-handle sites over the CORRECTED 103-handle incident pool.

Priced against 7 -> 12 (T's correction): the deliverable fails 7; with the 12 cofactors
zeroed it fails 12; the gap the cofactor freedom buys is 5.

Ordering: handles sorted by rt (how many target equations they reach), and 4-subsets taken
in lexicographic order over that ranking, so high-rt sites come first.  An interrupted run
is then a real measurement over a stated prefix.  C(103,4) = 4.4M is far beyond budget, so
a stopping point is expected and is reported.

Note on eq8680 / S = 0: no special handling is needed here.  The tuner measures the TRUE
score by re-propagation, so breaking S = 0 is counted as the collateral it is -- which is
exactly what my round-9 solve did by moving x_4432 and x_28730 by 3571 vs 3572 bits.
"""
import sys, os, json, time, itertools, collections, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import price as PR, fscore, sweep

BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 2400.0
D4 = [642, 28730, 29854, 31864]
prev = pickle.load(open('pricelead.pkl', 'rb'))
rt = {v['added']: v['rt'] for v in prev.values()}
for h in D4:
    rt.setdefault(h, 99)
POOL = sorted(rt, key=lambda u: (-rt[u], u))
print(f'pool {len(POOL)} handles, ranked by rt', flush=True)
print(f'  top 15: {[(u, rt[u]) for u in POOL[:15]]}', flush=True)

ctx = sweep.Ctx()
r = sweep.tune_fast(ctx, D4)
print(f'\ncalibration on the deliverable\'s four: {r["base_score"]} -> {r["score"]}  '
      f'{"PASSED" if r["score"] >= 39026 else "FAILED"}', flush=True)
if r['score'] < 39026:
    sys.exit(1)
print('reference points: deliverable 39026 (7 failing) ; cofactors zeroed 39021 (12 failing)',
      flush=True)

res = {}
best = (0, None)
t0 = time.time()
n = 0
stopped = None
for site in itertools.combinations(POOL, 4):
    n += 1
    try:
        rr = sweep.tune_fast(ctx, list(site), want=False)
    except Exception:
        continue
    if not rr.get('ok'):
        continue
    res[site] = (rr['score'], rr['nrows_target'])
    if rr['score'] > best[0]:
        best = (rr['score'], site)
        if rr['score'] > 39026:
            rr2 = sweep.tune_fast(ctx, list(site), want=True)
            if rr2.get('seed') is not None:
                v = rr2['eng'].forward(rr2['seed'])
                fn = f'M_e103_{rr2["score"]}.json'
                json.dump({f"x_{k}": int(v[k]) for k in range(PR.NV) if v[k] != 0},
                          open(fn, 'w'))
                print(f'  *** ABOVE BASELINE {rr["score"]} at {site} -> {fn} ***', flush=True)
    if n % 250 == 0:
        print(f'  [{n}] {time.time()-t0:.0f}s  best {best[0]}  last site {site}', flush=True)
        pickle.dump({'res': res, 'best': best, 'n': n}, open('enum103.pkl', 'wb'))
    if time.time() - t0 > BUDGET:
        stopped = site
        print(f'  [BUDGET reached after {n} sites; last site priced {site}]', flush=True)
        break

pickle.dump({'res': res, 'best': best, 'n': n, 'stopped_after': stopped},
            open('enum103.pkl', 'wb'))
sc = collections.Counter(v[0] for v in res.values())
print(f'\n=== DISTRIBUTION over {len(res)} priced 4-handle sites ===')
for k in sorted(sc, reverse=True):
    print(f'  {k}: {sc[k]}')
print(f'above 39026: {sum(1 for v in res.values() if v[0] > 39026)}')
print(f'equal 39026: {sum(1 for v in res.values() if v[0] == 39026)}')
print(f'at/below 39021 (the 12-failing far side): {sum(1 for v in res.values() if v[0] <= 39021)}')
print(f'priced out at 0 rows: {sum(1 for v in res.values() if v[1] == 0)}')
print(f'BEST {best[0]} {best[1]}')
print(f'STOPPING POINT: {n} of {len(list(itertools.combinations(range(len(POOL)), 4)))} '
      f'4-subsets, prefix in rt-ranked lexicographic order')
