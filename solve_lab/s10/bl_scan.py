"""bl_scan: generic cheap-prefilter scan over flip sets.

usage: bl_scan.py <frame:f2|canon> <tag> <mode:single|pair|triple> [prefix=2081,4287] [cands=cone|pins|all|<file>]
Writes bl_scan_<tag>.json = sorted [[score, nz, flips...]] and prints the top.
"""
import os, sys, json, time, itertools, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, cheap, FORBID

frame = sys.argv[1]; tag = sys.argv[2]; mode = sys.argv[3]
prefix = [int(x) for x in sys.argv[4].split(',') if x] if len(sys.argv) > 4 and sys.argv[4] != '-' else []
cs = sys.argv[5] if len(sys.argv) > 5 else 'cone'

if frame == 'f2':
    F = F2; v0 = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
else:
    F = CANON; v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
F.fwd(v0, rounds=6)
base, av0, nz0 = pot(v0)

RES = [a for a in range(L.NA) if av0[a]]
CONE = CANON.cone([21617, 29539]) | F.cone(RES)
BF = sorted((BOOL & CANON.FREE) - FORBID)
if cs == 'cone':
    cands = sorted(set(BF) & CONE)
elif cs == 'pins':
    LD = json.load(open(os.path.join(HERE, 'bl_pins3.json')))['loads']
    cands = sorted(set(b for _, b, _, _, _ in LD) & set(BF))
elif cs == 'all':
    cands = BF
elif cs == 'conepins':
    LD = json.load(open(os.path.join(HERE, 'bl_pins3.json')))['loads']
    cands = sorted((set(b for _, b, _, _, _ in LD) | (set(BF) & CONE)) & set(BF))
else:
    cands = sorted(json.load(open(os.path.join(HERE, cs))))

print(f'[{tag}] frame={frame} base={base[0]} nz={len(nz0)} {nz0[:8]}', flush=True)
print(f'[{tag}] prefix={prefix} cands={len(cands)} mode={mode}', flush=True)
if prefix:
    p, _, nz = cheap(v0, prefix, F)
    print(f'[{tag}] prefix alone: score {p[0]} nz {len(nz)} {nz[:10]}', flush=True)

if mode == 'single':
    sets = [(u,) for u in cands]
elif mode == 'pair':
    sets = list(itertools.combinations(cands, 2))
else:
    sets = list(itertools.combinations(cands, 3))
print(f'[{tag}] {len(sets)} combinations', flush=True)

t0 = time.time(); res = []; best = None
for i, s in enumerate(sets):
    p, v, nz = cheap(v0, list(prefix) + list(s), F)
    res.append((p[0], -p[1], list(s)))
    if best is None or p > best[0]:
        best = (p, list(s), v)
        if p[0] > 39026:
            T.save(v, os.path.join(HERE, f'bl_scan_{tag}_{p[0]}.json'))
            print(f'  *** {tag} {s}: score {p[0]} nz {nz}', flush=True)
    if i % 200 == 0:
        print(f'  {i}/{len(sets)} ({time.time()-t0:.0f}s) best {best[0][0]} via {best[1]}', flush=True)
res.sort(key=lambda r: (-r[0], r[1]))
json.dump([[int(a), int(b), [int(x) for x in c]] for a, b, c in res],
          open(os.path.join(HERE, f'bl_scan_{tag}.json'), 'w'))
print(f'\n[{tag}] TOP 30 (score, nz, flips):')
for r in res[:30]: print('   ', r)
print(f'[{tag}] best {best[0][0]} via {best[1]};  base {base[0]};  {time.time()-t0:.0f}s')
