"""bl_multi: beyond pairs/triples -- random k-subsets and a greedy beam over
boolean flip sets, cheap-filtered, then the engine on the shortlist."""
import os, sys, json, time, random, itertools
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, cheap, engine, FORBID
P = 2**256-2**32-977
random.seed(11)

frame = sys.argv[1] if len(sys.argv) > 1 else 'canon'
if frame == 'f2':
    F = F2; v0 = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
else:
    F = CANON; v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
F.fwd(v0, rounds=6)
base, av0, nz0 = pot(v0)
POOL = json.load(open(os.path.join(HERE, 'bl_cands78.json')))
PRIO = [438, 490, 5643, 5910, 6821, 12054, 12095, 13195, 14808, 16586, 16827,
        17406, 17760, 18022, 21074, 22562, 23751, 24365, 27393, 28005, 34974, 38625]
print(f'[multi/{frame}] base {base[0]} nz {len(nz0)} {nz0[:8]}  pool {len(POOL)} prio {len(PRIO)}', flush=True)

res = []
t0 = time.time()
# ---- random k-subsets ----
for pool, pn in ((POOL, 'cone76'), (PRIO, 'prio22')):
    for k in (2, 3, 5, 8, 13, 21):
        if k > len(pool): continue
        for _ in range(30):
            s = sorted(random.sample(pool, k))
            for pfn, pf in (('', []), ('b11', [4287])):
                p, v, nz = cheap(v0, pf + s, F)
                res.append((p[0], -p[1], pn, k, pfn, s))
        print(f'  {pn} k={k}: best so far {max(res)[0]} ({time.time()-t0:.0f}s)', flush=True)

# ---- greedy beam over the priority + cone pool ----
print('\n[greedy] beam over cone76 (width 3, depth 6)', flush=True)
beam = [((base[0], base[1]), [])]
seen = set()
for depth in range(6):
    cand = []
    for (sc, cur) in beam:
        for u in POOL:
            if u in cur: continue
            s = tuple(sorted(cur + [u]))
            if s in seen: continue
            seen.add(s)
            p, v, nz = cheap(v0, list(s), F)
            cand.append(((p[0], p[1]), list(s)))
    cand.sort(key=lambda t: -t[0][0] * 10**6 - t[0][1])
    beam = cand[:3]
    print(f'  depth {depth+1}: {[(c[0][0], -c[0][1], c[1]) for c in beam]} ({time.time()-t0:.0f}s)', flush=True)
    res += [(c[0][0], c[0][1], 'greedy', depth+1, '', c[1]) for c in beam]

res.sort(key=lambda r: (-r[0], r[1]))
json.dump([[int(a), int(b), c, int(d), e, [int(x) for x in f]] for a, b, c, d, e, f in res[:400]],
          open(os.path.join(HERE, f'bl_multi_{frame}.json'), 'w'))
print(f'\n[multi/{frame}] cheap top 20:')
for r in res[:20]: print('   ', r)

print(f'\n[multi/{frame}] engine on top 12', flush=True)
best = (base, None, v0)
for sc, nzn, pn, k, pfn, s in res[:12]:
    v = list(v0)
    for u in ([4287] if pfn == 'b11' else []) + s: v[u] = 1 - v[u]
    F.fwd(v, rounds=6)
    cur, vv, nz = engine(v, F, iters=50, budget=150)
    print(f'  {pn}/{pfn}/{s}: cheap {sc} -> ENGINE {cur[0]} nz {len(nz)} {nz[:8]}', flush=True)
    if cur > best[0]:
        best = (cur, (pn, pfn, s), vv)
        T.save(vv, os.path.join(HERE, f'bl_multi_{frame}_{cur[0]}.json'))
        if cur[0] > 39026:
            T.save(vv, os.path.join(HERE, 'bl_best.json')); print('  *** SAVED bl_best.json', flush=True)
print(f'[multi/{frame}] BEST {best[0][0]} via {best[1]}  (base {base[0]})  {time.time()-t0:.0f}s')
