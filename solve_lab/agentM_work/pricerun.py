"""Price L's candidate sites with handle values TUNED (not left at 0).

L's finding: the site fixes WHICH atoms break, the handle VALUES fix how many equations
they cost.  The deliverable's own site reads 13 failures with handles unset against a true
cost of 7, so an untuned score is meaningless.  Every candidate here is tuned before it is
scored, and the score is an exact re-propagation, never an estimate.
"""
import sys, os, json, time, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import price as PR

BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 150.0

# L's relayed candidates. Row 0 is the calibration row = the deliverable's own four.
CANDS = [
    ([642, 28730, 29854, 31864],   'c=27994 parent 4971.va   CALIBRATION (deliverable)', 13),
    ([10509, 20157, 32245, 33044], 'c=27634 parent 23762.va', 14),
    ([9541, 19546, 25227, 31891],  'c=31049 parent 14869.vb', 15),
    ([3260, 11588, 30400, 37248],  'c=30609 parent 14692.va', 15),
    ([2493, 3022, 6019, 15174],    'c=34711 parent 3154.vb',  15),
    ([1405, 3052, 4806, 16433],    'c=35056 parent 3154.va',  15),
    ([9337, 17894, 23336, 33996],  'c=6593  parent 25642.va', 15),
    ([19053, 21505, 22193, 23910], 'c=25642 parent 37991.vb', 15),
    ([10074, 16399, 16800, 35694], 'c=37034 parent 11182.va', 16),
    ([1768, 6389, 26662, 31362],   'c=23089 parent 37919.vb', 16),
    ([6254, 7439, 21115, 38560],   'c=37220 parent 33113.va', 16),
    ([1079, 15006, 15333, 32131],  'c=7383  parent 16012.vb', 16),
]

vd = PR.load_deliverable()
P = PR.TunedPricer(vd)

print('=== CALIBRATION, path 1: values supplied (the deliverable\'s own) ===', flush=True)
h0 = CANDS[0][0]
freed, demote = PR.closure(h0)
r = P.price_given(h0, {u: vd[u] for u in freed})
exp = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
nd = sum(1 for u in range(PR.NV) if r['v'][u] != vd[u])
print(f'  handles {h0} -> freed {freed} (5th DERIVED), demote {demote}')
print(f'  score {r["score"]}, fails {r["fails"]}, nbad {r["nbad"]}, vars differing {nd}')
print(f'  CALIBRATION(given) {"PASSED" if r["score"]==39026 and r["fails"]==exp and nd==0 else "FAILED"}',
      flush=True)

print('\n=== all candidates, path 2: handle values TUNED by greedy equation-space solve ===',
      flush=True)
results = []
t_all = time.time()
for hs, tag, incid in CANDS:
    t0 = time.time()
    try:
        rr = P.price_tuned(hs, budget=BUDGET, want=True)
    except Exception as ex:
        print(f'  {hs} {tag}: ERR {type(ex).__name__}: {ex}', flush=True)
        continue
    dt = time.time() - t0
    if not rr.get('ok'):
        print(f'  {hs} {tag}: {rr.get("why")}  ({dt:.0f}s)', flush=True)
        continue
    results.append((hs, tag, incid, rr))
    print(f'  {tag}\n     handles {hs} freed {len(rr["freed"])} knobs {rr["nknobs"]} '
          f'| L-incid {incid} | base(untuned) {rr["base_score"]} '
          f'-> TUNED {rr["score"]}  ({dt:.0f}s)', flush=True)
    if rr['score'] > 39026:
        eng = rr.get('eng'); sd = rr.get('seed')
        if eng is not None and sd is not None:
            v = eng.forward(sd)
            json.dump({f"x_{i}": int(v[i]) for i in range(PR.NV) if v[i] != 0},
                      open(f'M_site_{rr["score"]}.json', 'w'))
            print(f'     *** ABOVE BASELINE {rr["score"]} -> M_site_{rr["score"]}.json ***',
                  flush=True)

pickle.dump([(h, t, i, {k: v for k, v in r.items() if k not in ('eng', 'seed')})
             for h, t, i, r in results], open('pricerun.pkl', 'wb'))

print('\n=== SUMMARY (tuned) ===')
print(f'{"handles":34s} {"L-incid":>7s} {"untuned":>8s} {"TUNED":>7s}')
for hs, tag, incid, rr in sorted(results, key=lambda x: -x[3]['score']):
    print(f'{str(hs):34s} {incid:7d} {rr["base_score"]:8d} {rr["score"]:7d}')
tot = time.time() - t_all
if results:
    per = tot / len(results)
    print(f'\nthroughput: {per:.0f}s per candidate (tuned) -> ~{3600/per:,.0f}/hour single-core')
print('baseline to beat: 39026')
