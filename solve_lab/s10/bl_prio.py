"""bl_prio: pairs & triples over the PRIORITY boolean set = MUX controls that gate
a conditional constant LOAD whose pinned wire is an ancestor of a residual gadget
(a21617, a29539, a37662, a7930).  Cheap prefilter over every prefix branch, then
the full enriched engine on the shortlist.
"""
import os, sys, json, time, itertools, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, cheap, engine, FORBID
P = 2**256-2**32-977

PRIO = [438, 490, 5643, 5910, 6821, 12054, 12095, 13195, 14808, 16586, 16827,
        17406, 17760, 18022, 21074, 22562, 23751, 24365, 27393, 28005, 34974, 38625]
frame = sys.argv[1] if len(sys.argv) > 1 else 'canon'
if frame == 'f2':
    F = F2; v0 = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
else:
    F = CANON; v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
F.fwd(v0, rounds=6)
base, av0, nz0 = pot(v0)
print(f'[prio/{frame}] base {base[0]} nz {len(nz0)} {nz0[:10]}', flush=True)
print(f'[prio/{frame}] priority booleans: {PRIO} (all currently {sorted(set(v0[u] for u in PRIO))})', flush=True)

# triples restricted to the 12 gates whose pin lands directly in cone(a21617),
# cone(a29539) or cone(a7930)
TRI = [438, 5643, 6821, 12054, 13195, 16586, 16827, 17406, 21074, 22562, 23751, 24365]
PREFIX = {'': [], 'b00': [2081], 'b11': [4287], 'b01': [2081, 4287]}
combos = []
for pn, pf in PREFIX.items():
    for k in (1, 2):
        for s in itertools.combinations(PRIO, k):
            combos.append((pn, pf, list(s)))
for s in itertools.combinations(TRI, 3):
    combos.append(('', [], list(s)))
print(f'[prio/{frame}] {len(combos)} cheap evaluations', flush=True)

t0 = time.time(); res = []
for i, (pn, pf, s) in enumerate(combos):
    p, v, nz = cheap(v0, pf + s, F)
    res.append((p[0], -p[1], pn, s))
    if i % 200 == 0:
        print(f'  {i}/{len(combos)} ({time.time()-t0:.0f}s) best {max(res)[0]}', flush=True)
res.sort(key=lambda r: (-r[0], r[1]))
json.dump([[int(a), int(b), c, [int(x) for x in d]] for a, b, c, d in res],
          open(os.path.join(HERE, f'bl_prio_{frame}.json'), 'w'))
print(f'\n[prio/{frame}] cheap top 40:')
for r in res[:40]: print('   ', r)

# ---- engine on the shortlist ----
short, seen = [], set()
for sc, nznum, pn, s in res:
    key = (pn, tuple(s))
    if key in seen: continue
    seen.add(key); short.append((sc, nznum, pn, s))
    if len(short) >= 14: break
print(f'\n[prio/{frame}] running the enriched engine on {len(short)} shortlisted combos', flush=True)
best = (base, None, v0)
for sc, nznum, pn, s in short:
    v = list(v0)
    for u in PREFIX[pn] + s: v[u] = 1 - v[u]
    F.fwd(v, rounds=6)
    t1 = time.time()
    cur, vv, nz = engine(v, F, iters=50, budget=110)
    print(f'  {pn}+{s}: cheap {sc} -> ENGINE {cur[0]}  nz {len(nz)} {nz[:8]}  ({time.time()-t1:.0f}s)', flush=True)
    if cur > best[0]:
        best = (cur, (pn, s), vv)
        T.save(vv, os.path.join(HERE, f'bl_prio_{frame}_{cur[0]}.json'))
        print(f'   *** new best {cur[0]} via {pn}+{s}', flush=True)
        if cur[0] > 39026:
            T.save(vv, os.path.join(HERE, 'bl_best.json'))
            print('   *** SAVED bl_best.json', flush=True)
print(f'[prio/{frame}] BEST {best[0][0]} via {best[1]} (base {base[0]})  total {time.time()-t0:.0f}s')
