"""S12 step 14: is the COLLATERAL of a good activation repairable?

Apply the best-ranked activations, then run the repair engine from the activated
state and record the real integer score.  Activation is only worth anything if
the engine can win back more than the activation cost.
"""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
import engine as E
P = ad.P
BASE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE,'mod9118_0.json')
BUDGET = int(sys.argv[2]) if len(sys.argv) > 2 else 420
CANDS = sys.argv[3] if len(sys.argv) > 3 else 'ac_cands.json'
B = A.Base(BASE)
print(f'base {os.path.basename(BASE)} score {B.score0}', flush=True)
cands = json.load(open(os.path.join(HERE, CANDS)))
best = (B.score0, 'base')
for c in cands:
    act = {int(k): int(v_) for k, v_ in c['act'].items()}
    v = list(B.v0)
    for u, val in act.items(): v[u] = val
    if act: A.fwd_local(v, list(act))
    changed = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
    sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
    print(f'\n=== {c["tag"]}: activated score {sc} (lost {len(lost)}, gained {len(gained)}) ===', flush=True)
    try:
        v2, cur = E.run(v, f'ac_{c["tag"]}', iters=40, budget=BUDGET)
    except Exception as ex:
        print('  engine failed:', ex, flush=True); continue
    print(f'  {c["tag"]}: {sc} -> {cur[0]}', flush=True)
    if cur[0] > best[0]:
        best = (cur[0], c['tag'])
        T.save(v2, os.path.join(HERE, f'ac_best_{cur[0]}.json'))
    if cur[0] > 39026:
        T.save(v2, os.path.join(HERE, 'ac_best.json'))
        print('  *** BEAT 39026 ***', flush=True)
print(f'\nBEST after repair: {best}')
