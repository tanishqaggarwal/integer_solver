"""S10 step 35: the CORRECT global linear system.

A step delta must (a) zero the failing checks and (b) keep every SATISFIED check
satisfied.  Computing all 10,786 satisfied-check gradients is too slow, but we do
not need them: only checks actually REACHED by the step's support can change.

  1. support = union of the failing checks' gradient supports (non-boolean)
  2. for each u in support, forward-BFS the live (mod p) gate DAG -> reached checks
  3. build rows: failing checks -> -residual ; reached satisfied checks -> 0
  4. solve over GF(p) restricted to `support`

If this is consistent, the step fixes everything it can touch and breaks nothing.
If inconsistent, the certificate says exactly which combination is blocked.
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad
from newton import BOOL, solve_gfp

P = ad.P
atom_out = L.atom_out

v = L.load(os.path.join(HERE, 'forward_state.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
FAIL = [a for a in range(L.NA) if av[a] and a not in atom_out]
print('failing checks:', FAIL)


def live_reach(u):
    """Variables reachable from free input u through gates with nonzero mod-p
    sensitivity, plus the CHECK atoms thereby touched."""
    seen = {u}
    frontier = [u]
    checks = set()
    while frontier:
        nxt = []
        for w in frontier:
            for a in L.var_atoms[w]:
                if ad.dpart(a, w, vm) == 0:
                    continue
                if a not in atom_out:
                    checks.add(a)
                    continue
                t = atom_out[a][1]
                if t == w or t in seen:
                    continue
                if ad.dpart(a, t, vm) == 0:
                    continue
                seen.add(t); nxt.append(t)
        frontier = nxt
    return seen, checks


# 1. support
support = set()
grads = {}
for c in FAIL:
    g = {u: d for u, d in ad.grad(c, vm).items() if u not in BOOL}
    grads[c] = g
    support |= set(g)
support = sorted(support)
print(f'support: {len(support)} non-boolean free inputs -> {[f"x_{u}" for u in support]}')

# 2. reached checks
reached = set()
per_u = {}
for u in support:
    _, ch = live_reach(u)
    per_u[u] = ch
    reached |= ch
    print(f'  x_{u:<7} reaches {len(ch)} checks')
sat_reached = sorted(reached - set(FAIL))
print(f'\nreached checks: {len(reached)} total, {len(sat_reached)} currently satisfied')

# 3. rows
rows, rhs, names = [], [], []
for c in FAIL:
    rows.append(grads[c]); rhs.append((-av[c]) % P); names.append(('FAIL', c))
for c in sat_reached:
    g = {u: d for u, d in ad.grad(c, vm).items() if u in set(support)}
    if not g:
        continue
    rows.append(g); rhs.append(0); names.append(('SAT', c))
print(f'system: {len(rows)} rows x {len(support)} columns over GF(p)')

sol = solve_gfp(rows, rhs, support)
if sol is None:
    print('\n*** INCONSISTENT -- extracting the certificate')
    # find a minimal inconsistent subset by adding rows one at a time
    keep, krhs, kn = [], [], []
    for r, b, nm in zip(rows, rhs, names):
        keep.append(r); krhs.append(b); kn.append(nm)
        if solve_gfp(keep, krhs, support) is None:
            print(f'  becomes inconsistent when adding {nm}')
            print(f'  rows involved: {kn}')
            break
else:
    print(f'\n*** CONSISTENT: step over {len(sol)} inputs')
    for u, d in sol.items():
        v[u] = v[u] + d
    ad.fwd(v)
    av2 = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av2[a]]
    fail = L.failing_eqs(av2)
    print(f'   after step: nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}')
    T.save(v, os.path.join(HERE, 'constrained_step.json'))
