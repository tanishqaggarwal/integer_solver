"""bl_key4: the only booleans that can reach the failing equations.

The 24 atoms appearing in the 7 failing equations have a joint ancestor cone of
151 variables containing exactly TWO free booleans besides the MUX pair:
x_11368 and x_13195 (both are load-pin gates).  Enumerate all 16 assignments of
{x_2081, x_4287, x_11368, x_13195} in both frames, full engine on each.
"""
import os, sys, json, time, itertools
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, engine, FORBID
P = 2**256-2**32-977
KEY = [2081, 4287, 11368, 13195]
BUD = int(sys.argv[2]) if len(sys.argv) > 2 else 100
frame = sys.argv[1] if len(sys.argv) > 1 else 'f2'
if frame == 'f2':
    F = F2; v0 = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
else:
    F = CANON; v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
F.fwd(v0, rounds=6)
base, av0, nz0 = pot(v0)
print(f'[key4/{frame}] base {base[0]} nz {len(nz0)} {nz0[:8]}', flush=True)
print(f'[key4/{frame}] current key values: {[(u, v0[u]) for u in KEY]}', flush=True)

rows = []
best = (base, None, v0)
t0 = time.time()
for bits in itertools.product((0, 1), repeat=4):
    v = list(v0)
    for u, b in zip(KEY, bits): v[u] = b
    F.fwd(v, rounds=6)
    p, av, nz = pot(v)
    t1 = time.time()
    cur, vv, nz2 = engine(v, F, iters=60, budget=BUD)
    rows.append((cur[0], p[0], list(bits), len(nz2), nz2[:10]))
    print(f'  {bits}: plain {p[0]} (nz {len(nz)}) -> ENGINE {cur[0]} nz {len(nz2)} '
          f'{nz2[:10]}  ({time.time()-t1:.0f}s)', flush=True)
    if cur > best[0]:
        best = (cur, list(bits), vv)
        T.save(vv, os.path.join(HERE, f'bl_key4_{frame}_{cur[0]}.json'))
        print(f'   *** new best {cur[0]} at {bits}', flush=True)
        if cur[0] > 39026:
            T.save(vv, os.path.join(HERE, 'bl_best.json'))
            print('   *** SAVED bl_best.json', flush=True)
rows.sort(reverse=True)
json.dump(rows, open(os.path.join(HERE, f'bl_key4_{frame}.json'), 'w'))
print(f'[key4/{frame}] BEST {best[0][0]} at {best[1]} (base {base[0]})  {time.time()-t0:.0f}s', flush=True)
