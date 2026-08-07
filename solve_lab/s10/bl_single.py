"""bl_single: exhaustive single-flip cheap scan of ALL free booleans, in both frames.

Cheap = flip, forward-solve in the frame, count failing eqs + nonzero atoms.
This is the pre-filter feeding the pair/triple search.
"""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, cheap, FORBID

which = sys.argv[1] if len(sys.argv) > 1 else 'f2'
if which == 'f2':
    F = F2; v0 = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
else:
    F = CANON; v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
F.fwd(v0, rounds=6)
base, av0, nz0 = pot(v0)
print(f'[{which}] base {base[0]} nonzero {nz0}', flush=True)

BF = sorted((BOOL & F.FREE) - FORBID)
print(f'[{which}] free booleans: {len(BF)}', flush=True)

# cone membership (informational)
RES = [a for a in range(L.NA) if av0[a]]
CONE = F.cone(RES)
print(f'[{which}] residual cone vars {len(CONE)}, booleans in it {len(set(BF)&CONE)}', flush=True)

t0 = time.time(); res = []
for i, u in enumerate(BF):
    p, v, nz = cheap(v0, [u], F)
    res.append((p[0], -p[1], u, u in CONE))
    if i % 100 == 0:
        print(f'  {i}/{len(BF)} ({time.time()-t0:.0f}s) best {max(res)[0]}', flush=True)
res.sort(reverse=True)
json.dump([[int(a), int(b), int(c), bool(d)] for a, b, c, d in res],
          open(os.path.join(HERE, f'bl_single_{which}.json'), 'w'))
print(f'[{which}] base {base[0]}/{-base[1]}nz;  top 30 singles (score, nz, var, inCone):')
for r in res[:30]: print('   ', r)
print(f'[{which}] #flips scoring >= base: {sum(1 for r in res if r[0] >= base[0])}')
print(f'[{which}] #flips scoring > base : {sum(1 for r in res if r[0] > base[0])}')
print(f'total {time.time()-t0:.0f}s')
