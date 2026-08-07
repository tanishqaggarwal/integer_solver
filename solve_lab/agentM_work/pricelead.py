"""Price the lead: supersets of the deliverable's four, stage checks first.

Two corrections to my own earlier work are folded in here:

1. My "32 incident handles" was restricted to PRODUCT/BARE-defined variables -- the same
   blind spot I diagnosed for x_7068, now at scale.  30 of the truly incident handles are
   linearly defined and were invisible to it.  The pool below is over ALL definer forms.

2. There are two ways a handle reaches a target equation, and both count:
     (a) its DEFINER atom appears in a target equation   (corrupting it makes that atom nonzero)
     (b) some atom CONTAINING it appears in a target equation (its value moves that atom)
   The pool is the union.  (a)-only was 59; my old 32 was a shape-restricted (b).
"""
import sys, os, json, time, math, itertools, collections, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore
import sweep

D4 = [642, 28730, 29854, 31864]
STAGE = [23754, 35619, 9629]          # coordinator's priority order
T25 = json.load(open('baseline_sets.json'))['A']

# ---- atoms of the target equations, and rt per atom ----
rt_atom = collections.Counter()
for e in T25:
    for c, a in H.eqt[e][2]:
        if a >= 0 and c:
            rt_atom[a] += 1
A25 = set(rt_atom)

ATOM2VAR = {H.definer[u][0]: u for u in H.SEQ}
VAR2ATOM = {u: H.definer[u][0] for u in H.SEQ}

pool = {}
for u in H.SEQ:
    da = VAR2ATOM[u]
    via_a = rt_atom.get(da, 0)
    via_b = sum(rt_atom.get(a, 0) for a in H.occ.get(u, ()) if a in A25)
    if via_a or via_b:
        pool[u] = max(via_a, via_b)
print(f'CORRECTED POOL: {len(pool)} incident handles over all definer forms '
      f'(of {len(H.SEQ)} definer variables, {100.0*len(pool)/len(H.SEQ):.2f}%)')
old32 = set(json.load(open('incident_pool.json'))['incident_handles'])
print(f'  my earlier shape-restricted 32: {len(old32 & set(pool))} of 32 survive, '
      f'{len(set(pool) - old32)} handles were MISSING from it')
print(f'  deliverable four in pool: { {u: (u in pool) for u in D4} }')
print(f'  stage checks in pool    : { {u: (u in pool) for u in STAGE} }')

ctx = sweep.Ctx()

print('\n=== fast-tuner calibration (must reach 39,026 on the deliverable\'s four) ===',
      flush=True)
t0 = time.time()
r = sweep.tune_fast(ctx, D4, want=True)
print(f'  base {r["base_score"]} -> {r["score"]}   ({time.time()-t0:.1f}s, '
      f'{r["nrows_target"]} rows, {r["nknobs"]} knobs)', flush=True)
if r.get('seed') is not None:
    eng = r['eng']; v = eng.forward(r['seed'])
    full = fscore.score(eng.badatoms(v))
    print(f'  full re-propagation of same seed: {full}  incremental exact: {full==r["score"]}',
          flush=True)
if r['score'] < 39026:
    print('  CALIBRATION FAILED -- stopping, results would not be measurements.')
    sys.exit(1)
print('  CALIBRATION PASSED', flush=True)

rest = [u for u in pool if u not in D4]
# priority: the three stage checks first in the given order, then by rt
prio = STAGE + sorted((u for u in rest if u not in STAGE), key=lambda u: -pool[u])
sites = [tuple(sorted(D4 + [h])) for h in prio]

BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 2400.0
print(f'\n=== pricing {len(sites)} five-handle supersets (stage checks first) ===', flush=True)
res = {}
t_all = time.time()
best = (39026, None)
for i, s in enumerate(sites):
    added = [h for h in s if h not in D4][0]
    try:
        rr = sweep.tune_fast(ctx, list(s), want=True)
    except Exception as ex:
        print(f'  +x{added}: ERR {type(ex).__name__}', flush=True); continue
    if not rr.get('ok'):
        continue
    res[s] = {'added': added, 'rt': pool[added], 'score': rr['score'],
              'base': rr['base_score'], 'rows': rr['nrows_target'], 'knobs': rr['nknobs']}
    tag = '  <== STAGE CHECK' if added in STAGE else ''
    flag = '  *** ABOVE BASELINE ***' if rr['score'] > 39026 else ''
    print(f'  [{i+1}/{len(sites)}] +x{added:<6d} rt {pool[added]:2d} rows {rr["nrows_target"]:2d} '
          f'-> {rr["score"]}{tag}{flag}', flush=True)
    if rr['score'] > best[0]:
        best = (rr['score'], s)
        eng = rr['eng']; v = eng.forward(rr['seed'])
        fn = f'M_lead_{rr["score"]}.json'
        json.dump({f"x_{k}": int(v[k]) for k in range(PR.NV) if v[k] != 0}, open(fn, 'w'))
        print(f'     wrote {fn}', flush=True)
    if time.time() - t_all > BUDGET:
        print(f'  [budget reached after {i+1} sites]', flush=True); break

pickle.dump(res, open('pricelead.pkl', 'wb'))
sc = collections.Counter(v['score'] for v in res.values())
print(f'\n=== DISTRIBUTION over {len(res)} priced five-handle sites ===')
for k in sorted(sc, reverse=True):
    print(f'  {k}: {sc[k]}')
z = sum(1 for v in res.values() if v['rows'] == 0)
print(f'priced out at 0 rows once tuned: {z} of {len(res)}')
print(f'above 39026: {sum(1 for v in res.values() if v["score"] > 39026)}')
print(f'equal 39026: {sum(1 for v in res.values() if v["score"] == 39026)}')
print(f'BEST {best[0]}  {best[1]}')
