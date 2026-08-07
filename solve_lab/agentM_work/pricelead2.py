"""Six-handle sites: the deliverable's four plus every PAIR from the incident pool.

Five-handle supersets came back with nothing above 39,026 (89 of 98 exactly equal, 9 worse),
including all three stage checks.  Pairs are the next layer.  Pairs containing a stage check
are priced first so an interrupted run still answers the prioritised question.
"""
import sys, os, json, time, itertools, collections, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import price as PR, fscore, sweep

D4 = [642, 28730, 29854, 31864]
STAGE = [23754, 35619, 9629]
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 2400.0

prev = pickle.load(open('pricelead.pkl', 'rb'))
pool = sorted({v['added'] for v in prev.values()})
rt = {v['added']: v['rt'] for v in prev.values()}
print(f'pool {len(pool)} handles (excluding the deliverable\'s four)', flush=True)

ctx = sweep.Ctx()
r = sweep.tune_fast(ctx, D4)
print(f'calibration: {r["base_score"]} -> {r["score"]}  '
      f'{"PASSED" if r["score"]>=39026 else "FAILED"}', flush=True)
if r['score'] < 39026:
    sys.exit(1)

pairs = list(itertools.combinations(pool, 2))
st = set(STAGE)
pairs.sort(key=lambda p: (-len(st & set(p)), -(rt.get(p[0], 0) + rt.get(p[1], 0))))
print(f'{len(pairs)} pairs; {sum(1 for p in pairs if st & set(p))} contain a stage check',
      flush=True)

res = {}
best = (39026, None)
t0 = time.time()
for i, p in enumerate(pairs):
    site = tuple(sorted(D4 + list(p)))
    try:
        rr = sweep.tune_fast(ctx, list(site), want=(True))
    except Exception:
        continue
    if not rr.get('ok'):
        continue
    res[p] = {'score': rr['score'], 'rows': rr['nrows_target'], 'knobs': rr['nknobs'],
              'base': rr['base_score']}
    if rr['score'] > best[0]:
        best = (rr['score'], p)
        eng = rr['eng']; v = eng.forward(rr['seed'])
        fn = f'M_lead6_{rr["score"]}.json'
        json.dump({f"x_{k}": int(v[k]) for k in range(PR.NV) if v[k] != 0}, open(fn, 'w'))
        print(f'  *** ABOVE BASELINE {rr["score"]} at +{p} -> {fn} ***', flush=True)
    if i % 200 == 0:
        print(f'  [{i}/{len(pairs)}] {time.time()-t0:.0f}s  best {best[0]}', flush=True)
    if time.time() - t0 > BUDGET:
        print(f'  [budget after {i+1} pairs]', flush=True); break

pickle.dump(res, open('pricelead2.pkl', 'wb'))
sc = collections.Counter(v['score'] for v in res.values())
print(f'\n=== DISTRIBUTION over {len(res)} priced six-handle sites ===')
for k in sorted(sc, reverse=True):
    print(f'  {k}: {sc[k]}')
print(f'priced out at 0 rows: {sum(1 for v in res.values() if v["rows"]==0)} of {len(res)}')
print(f'above 39026: {sum(1 for v in res.values() if v["score"]>39026)}')
print(f'equal 39026: {sum(1 for v in res.values() if v["score"]==39026)}')
print(f'BEST {best[0]} {best[1]}')
