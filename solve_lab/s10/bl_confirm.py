"""bl_confirm: engine on the pair78-specific top pairs not already engined.
Chunked + JSONL-checkpointed so a kill costs at most one candidate."""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import CANON, F2, pot, engine
OUT = os.path.join(HERE, 'bl_confirm.jsonl')
CAND = [[2455, 33287], [24517, 33287], [11368, 22562], [11368, 24365],
        [23751, 24365], [438, 16827], [12054, 16586, 24365], [16586, 24365, 17406]]
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        try: done.add(tuple(json.loads(line)['flips']))
        except Exception: pass
v0 = L.load(os.path.join(HERE, 'mod9118_0.json')); CANON.fwd(v0, rounds=6)
base, _, nz0 = pot(v0)
print(f'base {base[0]} nz {nz0}', flush=True)
for s in CAND:
    if tuple(s) in done:
        print(f'  skip {s} (checkpointed)', flush=True); continue
    v = list(v0)
    for u in s: v[u] = 1 - v[u]
    CANON.fwd(v, rounds=6)
    p0, _, _ = pot(v)
    t0 = time.time()
    cur, vv, nz = engine(v, CANON, iters=50, budget=70)
    rec = {'flips': s, 'cheap': p0[0], 'engine': cur[0], 'nz': nz, 'secs': round(time.time()-t0)}
    with open(OUT, 'a') as f: f.write(json.dumps(rec) + '\n'); f.flush()
    print(f'  {s}: cheap {p0[0]} -> ENGINE {cur[0]} nz {len(nz)} {nz[:8]}', flush=True)
    if cur[0] > 39026:
        T.save(vv, os.path.join(HERE, 'bl_best.json')); print('  *** SAVED bl_best.json', flush=True)
print('done', flush=True)
