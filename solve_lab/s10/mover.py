"""S11 step 95: is (A, B) reachable by ANY boolean flip?

neutral.py settled the score-neutral bits: all 300 tested are completely inert -- they
do not move A or B at all, so the free branch directions never reach the residual.
The complementary question is the decisive one: does ANY boolean free input move A or
B, at any price?  If none does, then on this branch A and B are functions of the seven
literal constants alone and no boolean setting can change them.

Usage: mover.py [state.json] [N]
"""
import os, sys, json, time, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
from chunk import load
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PIN_39013.json'
N = int(sys.argv[2]) if len(sys.argv) > 2 else 500
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
BASE = L.NEQ - len(L.failing_eqs(L.all_atom_values(base)))
A0, B0 = base[35389] % P, base[6671] % P
rs = load('bool_' + os.path.basename(src).replace('.json', ''))
COST = sorted([r['u'] for r in rs if r.get('score', -1) < BASE])
print('%s: %d census results, %d bits that COST score' % (src, len(rs), len(COST)),
      flush=True)
moves, inert, t0 = [], 0, time.time()
sc = collections.Counter()
for i, u in enumerate(COST[:N]):
    v = list(base)
    v[u] = 1 - v[u]
    ad.fwd(v, rounds=6)
    dA, dB = (v[35389] % P - A0) % P, (v[6671] % P - B0) % P
    s = L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))
    if dA or dB:
        moves.append((u, s, str(dA), str(dB)))
        sc[s] += 1
    else:
        inert += 1
    if i % 50 == 0:
        print('  %d/%d (%.0fs) moving %d inert %d'
              % (i, min(N, len(COST)), time.time() - t0, len(moves), inert), flush=True)
print('\ncosting bits that MOVE (A, B): %d;  inert: %d' % (len(moves), inert))
if moves:
    moves.sort(key=lambda m: -m[1])
    print('best-scoring movers:')
    for u, s, a, b in moves[:12]:
        print('   x%-6d score %-6d dA %s...  dB %s...' % (u, s, a[:24], b[:24]))
json.dump({'A0': str(A0), 'B0': str(B0), 'moves': moves},
          open(os.path.join(HERE, 'movers.json'), 'w'))
print('saved movers.json')
