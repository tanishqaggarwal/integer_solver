"""bl_final: run the full enriched engine on the global shortlist gathered from
every cheap scan (bl_scan_*.json, bl_prio_*.json, bl_multi_*.json)."""
import os, sys, json, glob, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, engine, FORBID

PREFIX = {'c_b11_all': [4287], 'c_b00': [2081], 'c_b01': [2081, 4287],
          'f2_b11': [4287], 'f2_b00': [2081], 'f2_b01': [2081, 4287],
          'c_pair78': [], 'c_pair_cone': [], 'f2_pair_pins': []}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
BUD = int(sys.argv[2]) if len(sys.argv) > 2 else 200

cands = []
for f in sorted(glob.glob(os.path.join(HERE, 'bl_scan_*.json'))):
    tag = os.path.basename(f)[8:-5]
    if tag not in PREFIX: continue
    frame = 'f2' if tag.startswith('f2') else 'canon'
    rows = json.load(open(f))
    for sc, nz, flips in rows[:14]:
        cands.append((sc, nz, frame, PREFIX[tag] + list(flips), tag))
for f in sorted(glob.glob(os.path.join(HERE, 'bl_prio_*.json'))):
    frame = 'f2' if 'f2' in os.path.basename(f) else 'canon'
    PF = {'': [], 'b00': [2081], 'b11': [4287], 'b01': [2081, 4287]}
    for sc, nz, pn, flips in json.load(open(f))[:14]:
        cands.append((sc, nz, frame, PF[pn] + list(flips), os.path.basename(f)))
cands.sort(key=lambda r: (-r[0], r[1]))
seen, short = set(), []
for c in cands:
    k = (c[2], tuple(sorted(c[3])))
    if k in seen: continue
    seen.add(k); short.append(c)
    if len(short) >= N: break
print(f'shortlist {len(short)} of {len(cands)} scanned candidates', flush=True)

V = {}
V['canon'] = L.load(os.path.join(HERE, 'mod9118_0.json')); CANON.fwd(V['canon'], rounds=6)
V['f2'] = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json')); F2.fwd(V['f2'], rounds=6)
BASE = {k: pot(v)[0] for k, v in V.items()}
print(f'bases: {BASE}', flush=True)
best = (max(BASE.values()), None)
for sc, nz, frame, flips, tag in short:
    F = F2 if frame == 'f2' else CANON
    v = list(V[frame])
    for u in flips: v[u] = 1 - v[u]
    F.fwd(v, rounds=6)
    t0 = time.time()
    cur, vv, nzz = engine(v, F, iters=60, budget=BUD)
    print(f'  [{frame}] {tag} {flips}: cheap {sc} -> ENGINE {cur[0]} nz {len(nzz)} '
          f'{nzz[:8]} ({time.time()-t0:.0f}s)', flush=True)
    if cur[0] > best[0]:
        best = (cur[0], (frame, flips))
        T.save(vv, os.path.join(HERE, f'bl_final_{cur[0]}.json'))
        if cur[0] > 39026:
            T.save(vv, os.path.join(HERE, 'bl_best.json'))
            print('  *** SAVED bl_best.json', flush=True)
print(f'FINAL BEST {best[0]} via {best[1]}', flush=True)
