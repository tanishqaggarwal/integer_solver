"""Enumerate 4-handle sites with the calibrated incremental engine.

rt-ranked lexicographic order, so an interrupted run is a real measurement over a stated
prefix.  Checkpoints frequently and is safe to kill at any moment -- the campaign may need
this agent as a verifier at short notice.
"""
import sys, os, json, time, itertools, collections, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore, price as PR

BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 1500.0
D4 = [642, 28730, 29854, 31864]
prev = pickle.load(open('pricelead.pkl', 'rb'))
rt = {v['added']: v['rt'] for v in prev.values()}
for h in D4:
    rt.setdefault(h, 99)
POOL = sorted(rt, key=lambda u: (-rt[u], u))
TOTAL = len(list(itertools.combinations(range(len(POOL)), 4)))
print(f'pool {len(POOL)}, total 4-subsets {TOTAL:,}', flush=True)

cal = ieng.tune(D4)
print(f'calibration: {cal["base_score"]} -> {cal["score"]}  '
      f'{"PASSED" if cal["score"] >= 39026 else "FAILED"}', flush=True)
if cal['score'] < 39026:
    sys.exit(1)
print('reference: deliverable 39026 (7 failing) ; cofactors zeroed 39021 (12 failing)',
      flush=True)

dist = collections.Counter()
best = (0, None)
above = []
n = 0
t0 = time.time()
last = None
for s in itertools.combinations(POOL, 4):
    n += 1
    try:
        r = ieng.tune(list(s))
    except Exception:
        continue
    if not r.get('ok'):
        continue
    sc = r['score']
    dist[sc] += 1
    last = s
    if sc > best[0]:
        best = (sc, s)
    if sc > 39026:
        above.append((s, sc))
        r2 = ieng.tune(list(s), want=True)
        if r2.get('changes'):
            bad, v = ieng.resid(ieng.V_UNC, ieng.BAD_UNC, r2['changes'], r2['pin'])
            fn = f'M_isweep_{sc}.json'
            json.dump({f"x_{k}": int(v[k]) for k in range(ieng.NV) if v[k] != 0}, open(fn, 'w'))
            print(f'  *** ABOVE BASELINE {sc} at {s} -> {fn} ***', flush=True)
    if n % 5000 == 0:
        el = time.time() - t0
        print(f'  [{n:,}/{TOTAL:,}] {el:.0f}s  {n/el:.0f} sites/s  best {best[0]}', flush=True)
        pickle.dump({'dist': dict(dist), 'best': best, 'n': n, 'last': last,
                     'above': above}, open('isweep.pkl', 'wb'))
    if time.time() - t0 > BUDGET:
        print(f'  [budget after {n:,} sites]', flush=True)
        break

el = time.time() - t0
pickle.dump({'dist': dict(dist), 'best': best, 'n': n, 'last': last, 'above': above},
            open('isweep.pkl', 'wb'))
print(f'\n=== DISTRIBUTION over {n:,} priced 4-handle sites ({el:.0f}s, '
      f'{n/el:.0f} sites/s) ===')
for k in sorted(dist, reverse=True):
    print(f'  {k}: {dist[k]:,}')
print(f'above 39026: {sum(v for k, v in dist.items() if k > 39026)}')
print(f'equal 39026: {dist.get(39026, 0):,}')
print(f'below 39021: {sum(v for k, v in dist.items() if k < 39021):,}')
print(f'BEST {best[0]} at {best[1]}')
print(f'STOPPING POINT: {n:,} of {TOTAL:,} in rt-ranked lexicographic order; '
      f'last site priced {last}')
