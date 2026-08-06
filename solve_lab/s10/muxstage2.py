"""S10 step 10: greedy equation-space repair on the MUX branch."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tools as T
import lib as L

v = T.L.load(os.path.join(HERE, 'mux_on.json'))
av, nz, fail = T.measure(v, 'start')

best = (L.NEQ - len(fail), list(v))
for rnd in range(14):
    av, nz, fail = T.measure(v, f'round {rnd}')
    if not nz:
        print('ALL ATOMS ZERO'); break
    # pick the move that maximises the equation-space score across all nonzero atoms
    cands = []
    for a in nz:
        for sc, u, nvl in T.try_fix(v, a)[:3]:
            cands.append((sc, a, u, nvl))
    if not cands:
        print('  no move available'); break
    cands.sort(key=lambda t: -t[0])
    sc, a, u, nvl = cands[0]
    cur = L.NEQ - len(fail)
    print(f'  best move: fix a{a} via x_{u} -> score {sc} (was {cur})')
    if sc <= cur and rnd > 0:
        print('  no improving move; stopping'); break
    L.ripple(v, {u: nvl})
    if sc > best[0]:
        best = (sc, list(v))

print(f'\nbest score reached: {best[0]}')
T.save(best[1], os.path.join(HERE, 'mux_greedy.json'))
