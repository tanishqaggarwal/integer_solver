"""S11 step 94: do the SCORE-NEUTRAL boolean flips move A and B?

boolcensus finds many boolean free inputs whose flip costs nothing -- the score stays
at the 39,013 attractor and the residual stays the same five checks.  Individually
none of them zeroes A or B, but "does not zero it" and "does not move it" are very
different facts, and only the second is a dead end.  So record the actual values.

If a neutral bit moves (A, B) at all, then the achievable set under flips of the
neutral bits is a sum of contributions and the question becomes a two-target subset
problem mod p -- a real target.  If every neutral bit leaves (A, B) untouched, the
neutral directions are genuinely inert and the branch is not reachable that way.

Usage: neutral.py [state.json] [MAX]
"""
import os, sys, json, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
from chunk import load
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PIN_39013.json'
MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 400
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
BASE = L.NEQ - len(L.failing_eqs(L.all_atom_values(base)))
A0, B0 = base[35389] % P, base[6671] % P
rs = load('bool_' + os.path.basename(src).replace('.json', ''))
NEU = [r['u'] for r in rs if r.get('score') == BASE]
print('%s: score %d; %d census results, %d score-neutral bits'
      % (src, BASE, len(rs), len(NEU)), flush=True)
print('A0 = %d\nB0 = %d' % (A0, B0), flush=True)

moves, inert, t0 = [], 0, time.time()
for i, u in enumerate(NEU[:MAX]):
    v = list(base)
    v[u] = 1 - v[u]
    ad.fwd(v, rounds=6)
    dA, dB = (v[35389] % P - A0) % P, (v[6671] % P - B0) % P
    if dA or dB:
        moves.append((u, dA, dB))
    else:
        inert += 1
    if i % 50 == 0:
        print('  %d/%d  (%.0fs)  moving %d, inert %d'
              % (i, min(MAX, len(NEU)), time.time() - t0, len(moves), inert),
              flush=True)
print('\nscore-neutral bits that MOVE (A, B): %d;  inert: %d'
      % (len(moves), inert), flush=True)
for u, dA, dB in moves[:20]:
    print('   x%-6d dA = %s...  dB = %s...' % (u, str(dA)[:28], str(dB)[:28]))
json.dump({'A0': str(A0), 'B0': str(B0),
           'moves': [[u, str(a), str(b)] for u, a, b in moves]},
          open(os.path.join(HERE, 'neutral.json'), 'w'))
print('saved neutral.json')

if moves:
    # is the contribution ADDITIVE?  test a few pairs against the sum of singles
    print('\nadditivity test on pairs:', flush=True)
    for (u1, a1, b1), (u2, a2, b2) in zip(moves[:4], moves[1:5]):
        v = list(base)
        v[u1] = 1 - v[u1]
        v[u2] = 1 - v[u2]
        ad.fwd(v, rounds=6)
        dA, dB = (v[35389] % P - A0) % P, (v[6671] % P - B0) % P
        print('   x%d+x%d: additive %s' % (u1, u2,
              (dA == (a1 + a2) % P and dB == (b1 + b2) % P)), flush=True)
