"""Exhaustive enumeration over subsets of the INCIDENT p-handle set.

A subset S = the set of p-handle atoms allowed to be nonzero, i.e. the handles whose
defining relation h = p*u is broken.  The 39,026 witness is exactly the subset
{642, 28730, 29854, 31864} of size 4, so the enumeration contains the known answer and the
calibration is built in.

Run in increasing order of space so a partial run is still a result:
    2^12 = 4,096    incident against the deliverable's own 7 failures
    2^16 = 65,536   incident against T's 12-equation far side  (what I was told to price against)
    2^18 = 262,144  incident against my 25-equation uncorrupted baseline

No ranking, no truncation, no early cutoff: per Q, an atom can be nonzero inside an equation
that still sums to zero, so incidence filters REACHABILITY, not cost, and a subset's price
cannot be bounded below by its incidence.  Every subset is priced by re-propagation.
"""
import sys, os, json, time, itertools, collections, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng, fscore, price as PR

WHICH = sys.argv[1] if len(sys.argv) > 1 else '12'
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 6000.0

PF = json.load(open('pfamily.json'))
H12 = sorted({v['h'] for v in PF['incident_7'].values()})
H16 = sorted({v['h'] for v in PF['incident_12'].values()})
H18 = sorted({v['h'] for v in PF['incident_25'].values()})
SETS = {'12': H12, '16': H16, '18': H18}
HL = SETS[WHICH]
D4 = [642, 28730, 29854, 31864]

print(f'incident handle sets: |H12|={len(H12)} |H16|={len(H16)} |H18|={len(H18)}', flush=True)
print(f'enumerating subsets of H{WHICH}: {HL}', flush=True)
print(f'  deliverable {D4} inside this set: {set(D4) <= set(HL)}', flush=True)
print(f'  total subsets 2^{len(HL)} = {2**len(HL):,}', flush=True)

cal = ieng.tune(D4)
print(f'calibration on the witness subset: {cal["base_score"]} -> {cal["score"]}  '
      f'{"PASSED" if cal["score"] >= 39026 else "FAILED"}', flush=True)
if cal['score'] < 39026:
    sys.exit(1)

dist = collections.Counter()
bysize = collections.defaultdict(collections.Counter)
best = (0, None)
above = []
n = 0
t0 = time.time()
last = None
# increasing |S| so a partial run is a clean statement about small supports
order = [S for k in range(len(HL) + 1) for S in itertools.combinations(HL, k)]
for S in order:
    n += 1
    try:
        r = ieng.tune(list(S)) if S else {'ok': True, 'score': ieng.NEQ - len(ieng.FAILS_UNC)}
    except Exception:
        continue
    if not r.get('ok'):
        continue
    sc = r['score']
    dist[sc] += 1
    bysize[len(S)][sc] += 1
    last = S
    if sc > best[0]:
        best = (sc, S)
    if sc > 39026:
        above.append((S, sc))
        r2 = ieng.tune(list(S), want=True)
        if r2.get('changes'):
            bad, v = ieng.resid(ieng.V_UNC, ieng.BAD_UNC, r2['changes'], r2['pin'])
            fn = f'M_sub{WHICH}_{sc}_{"_".join(map(str, S))}.json'[:120]
            json.dump({f"x_{k}": int(v[k]) for k in range(ieng.NV) if v[k] != 0}, open(fn, 'w'))
            print(f'  *** ABOVE BASELINE {sc} at S={S} -> {fn} ***', flush=True)
    if n % 2000 == 0:
        el = time.time() - t0
        print(f'  [{n:,}/{2**len(HL):,}] {el:.0f}s {n/el:.0f}/s  best {best[0]}  |S|={len(S)}',
              flush=True)
        pickle.dump({'dist': dict(dist), 'bysize': {k: dict(v) for k, v in bysize.items()},
                     'best': best, 'n': n, 'last': last, 'above': above},
                    open(f'enumsub{WHICH}.pkl', 'wb'))
    if time.time() - t0 > BUDGET:
        print(f'  [budget after {n:,}]', flush=True); break

el = time.time() - t0
pickle.dump({'dist': dict(dist), 'bysize': {k: dict(v) for k, v in bysize.items()},
             'best': best, 'n': n, 'last': last, 'above': above, 'complete': n == 2**len(HL)},
            open(f'enumsub{WHICH}.pkl', 'wb'))
print(f'\n=== DISTRIBUTION over {n:,} subsets of H{WHICH} '
      f'({el:.0f}s, {n/max(el,1e-9):.0f}/s, complete={n == 2**len(HL)}) ===')
for k in sorted(dist, reverse=True):
    print(f'  {k}: {dist[k]:,}')
print(f'\nabove 39026: {sum(v for k, v in dist.items() if k > 39026)}')
print(f'equal 39026: {dist.get(39026, 0):,}')
print(f'BEST {best[0]} at S={best[1]}')
print(f'\nbest score by support size |S|:')
for k in sorted(bysize):
    mx = max(bysize[k])
    print(f'  |S|={k:2d}: {sum(bysize[k].values()):7,} subsets, best {mx}, '
          f'count at best {bysize[k][mx]:,}')
