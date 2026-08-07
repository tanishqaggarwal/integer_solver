"""The only cheap single-atom growth: absorb a23618 (frees x_28730, adds exactly 1 equation)."""
import sys, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/grow23618.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)
    sys.stdout.flush()


R = G.R0 + [23618]
P = G.private_vars(R)
Eqs = sorted({e for x in R for e in G.EQCO[x]})
say('R =', sorted(R))
say('private vars =', P)
say('|E| =', len(Eqs), Eqs)
t0 = time.time()
r = G.evaluate(R, verbose=False)
if r is None:
    say('NONLINEAR')
    raise SystemExit
cost, ne, k, S, sol, PP, RR = r
say('maxsat = %d of %d  -> COST %d  (%.0fs)' % (k, ne, cost, time.time() - t0))
say('satisfied subset:', S)
w = G.realise(RR, PP, sol)
ff = E.eqfails(E.badatoms(w))
say('EXACT: %d failing equations -> score %d' % (len(ff), 39033 - len(ff)))
say('failing:', sorted(ff))
out = '%s/grow23618_%d.json' % (OD, 39033 - len(ff))
json.dump({f"x_{i}": str(int(w[i])) for i in range(E.NV) if w[i] != 0}, open(out, 'w'))
say('wrote', out)
for u in PP:
    say('  x_%d: %db -> %db' % (u, abs(G.V0[u]).bit_length(), abs(w[u]).bit_length()))
