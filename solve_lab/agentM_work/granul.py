"""Rule 9 applied to my own enumeration: vary the solver's granularity.

The enumeration is exhaustive over SUBSETS, but each subset's price is the output of ONE
solver at ONE setting -- `ieng.tune` with nprobe=10, budget=30s, which greedily accumulates
rows and then probes 10 points along the solution chain.  So "nothing above 39,026" is, per
rule 9, a statement about that granularity until the granularity is varied.

Test: re-price a sample at nprobe=80 / budget=180 (8x the probes, 6x the time) and compare
score-for-score against the default.  Sample is deliberately two-part:

  A. every superset of the witness with |W| <= 8 inside H16 -- the region that actually
     contains all 114 optima, so if extra probing buys anything it should show here;
  B. a seeded uniform random sample of the whole 2^16, so the test is not confined to the
     region I already believe is good.

Any subset whose score MOVES is the finding; a subset that moves ABOVE 39,026 is the win.
"""
import sys, os, json, time, itertools, random, collections

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR); sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
import ieng                                                    # noqa: E402

W4 = (642, 28730, 29854, 31864)
PF = json.load(open('pfamily.json'))
H16 = sorted({v['h'] for v in PF['incident_12'].values()})
REST = [h for h in H16 if h not in W4]

A = [tuple(sorted(W4 + c)) for k in range(0, 5)
     for c in itertools.combinations(REST, k)]
rng = random.Random(9)
ALL = [W for k in range(len(H16) + 1) for W in itertools.combinations(H16, k)]
B = rng.sample(ALL, 1000)
SAMPLE = list(dict.fromkeys(A + B))
print(f'sample: {len(A)} witness-supersets (|W|<=8) + 1000 uniform = {len(SAMPLE)} distinct',
      flush=True)

moved = []; up = []
t0 = time.time()
delta = collections.Counter()
for n, W in enumerate(SAMPLE, 1):
    lo = ieng.tune(list(W))
    hi = ieng.tune(list(W), nprobe=80, budget=180.0)
    if not (lo.get('ok') and hi.get('ok')):
        continue
    d = hi['score'] - lo['score']
    delta[d] += 1
    if d:
        moved.append((W, lo['score'], hi['score']))
    if hi['score'] > 39026:
        up.append((W, hi['score']))
        r2 = ieng.tune(list(W), nprobe=80, budget=180.0, want=True)
        if r2.get('changes'):
            bad, v = ieng.resid(ieng.V_UNC, ieng.BAD_UNC, r2['changes'], r2['pin'])
            fn = f'M_gran_{hi["score"]}_{"_".join(map(str, W))}.json'[:120]
            json.dump({f"x_{k}": int(v[k]) for k in range(ieng.NV) if v[k] != 0}, open(fn, 'w'))
            print(f'  *** ABOVE 39026 at HIGH GRANULARITY: {hi["score"]} W={W} -> {fn} ***',
                  flush=True)
    if n % 200 == 0:
        print(f'  [{n}/{len(SAMPLE)}] {time.time()-t0:.0f}s  moved {len(moved)}  above {len(up)}',
              flush=True)

print(f'\n=== granularity test, {len(SAMPLE)} subsets, {time.time()-t0:.0f}s ===')
print(f'score(nprobe=80,budget=180) - score(nprobe=10,budget=30):')
for d in sorted(delta):
    print(f'  delta {d:+d}: {delta[d]:,}')
print(f'subsets whose score MOVED : {len(moved)}')
print(f'subsets now above 39,026  : {len(up)}')
for W, a, b in moved[:20]:
    print(f'    {W}  {a} -> {b}')
json.dump({'n': len(SAMPLE), 'delta_hist': {str(k): v for k, v in delta.items()},
           'moved': [[list(w), a, b] for w, a, b in moved],
           'above': [[list(w), s] for w, s in up]},
          open('granul.json', 'w'), indent=1)
