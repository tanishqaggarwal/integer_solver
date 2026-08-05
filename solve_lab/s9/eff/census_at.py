"""Effective-footprint knob census at an arbitrary baseline.
Usage: python3 census_at.py <assignment.json> <tag>"""
import pickle, sys, time, collections
import lib as L, model as MD

path = sys.argv[1]
tag = sys.argv[2]
v0 = L.load(path)
MD.BASEP = [MD.prim_val(a, v0) for a in range(L.NA)]
base_av = L.all_atom_values(v0)
base_eq = [L.eq_value(i, base_av) for i in range(L.NEQ)]
fails0 = [i for i, x in enumerate(base_eq) if x]
print(f'{tag}: {L.NEQ-len(fails0)}/{L.NEQ}  failing={fails0}')

DELTA = 1234567891011
foot_eq = {}
foot_at = {}
rev = collections.defaultdict(set)
t0 = time.time()
for x in range(L.NVARS):
    v = list(v0)
    ch, st = L.ripple(v, {x: v0[x] + DELTA})
    cand = set()
    for u in ch:
        cand.update(MD.prim_var_atoms[u])
    tou = {a: MD.prim_val(a, v) for a in cand}
    tou = frozenset(a for a in tou if tou[a] != MD.BASEP[a])
    foot_at[x] = tou
    for a in tou:
        rev[a].add(x)
    # equation footprint: recompute only affected equations, exactly
    av = list(base_av)
    for u in ch:
        for a in L.var_atoms[u]:
            av[a] = L.evalpoly(L.polys[a], v)
    eqs = frozenset(i for i in L.eqs_of_atoms(set(tou) | set(a for u in ch for a in L.var_atoms[u]))
                    if L.eq_value(i, av) != base_eq[i])
    foot_eq[x] = eqs
print(f'census done {time.time()-t0:.0f}s')

gain = sum(1 for x in range(L.NVARS) if len(foot_eq[x]) < len(L.var_eqs.get(x, ())))
removed = sum(len(L.var_eqs.get(x, ())) - len(foot_eq[x]) for x in range(L.NVARS)
              if len(foot_eq[x]) < len(L.var_eqs.get(x, ())))
noop = sum(1 for x in range(L.NVARS) if not foot_eq[x])
syn_tot = sum(len(L.var_eqs.get(x, ())) for x in range(L.NVARS))
eff_tot = sum(len(foot_eq[x]) for x in range(L.NVARS))
print(f'  syntactic footprint total = {syn_tot}')
print(f'  effective  footprint total = {eff_tot}   ({100.0*(syn_tot-eff_tot)/syn_tot:.1f}% smaller)')
print(f'  variables that STRICTLY gain = {gain}   equations removed = {removed}')
print(f'  variables that are pure no-ops (effective footprint empty) = {noop}')
h = collections.Counter(len(foot_eq[x]) for x in range(L.NVARS))
print('  size histogram <=12:', {k: h[k] for k in sorted(h) if k <= 12})
pickle.dump({'foot_eq': foot_eq, 'foot_at': foot_at, 'rev': dict(rev), 'base': MD.BASEP,
             'fails0': fails0}, open(f'census_{tag}.pkl', 'wb'))
print(f'wrote census_{tag}.pkl')
