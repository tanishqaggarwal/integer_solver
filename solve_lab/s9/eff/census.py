"""Effective-footprint knob census over ALL 38,748 variables, at a given baseline.

For each variable x we perturb x by a large random delta, run the gate-preserving
forward ripple, and record EXACTLY which atoms change value and hence which
equations change value.  That is the *effective* footprint, versus the *syntactic*
footprint (equations whose text mentions x).

Usage:  python3 census.py [24|22]
Writes: census24.pkl / census22.pkl
"""
import pickle, random, sys, time, collections
import lib as L

which = sys.argv[1] if len(sys.argv) > 1 else '24'
base_path = L.BEST24 if which == '24' else L.BEST22
v0 = L.load(base_path)
av0 = L.all_atom_values(v0)
base_eq = {}
for i in range(L.NEQ):
    base_eq[i] = L.eq_value(i, av0)
fails0 = [i for i, x in base_eq.items() if x]
print(f'baseline {base_path}: {L.NEQ - len(fails0)}/{L.NEQ}, failing={fails0}')

DELTAS = [1234567891011, 98765432109876543, 7]
random.seed(12345)

res = {}
t0 = time.time()
for x in range(L.NVARS):
    eq_union = set()
    at_union = set()
    per_delta = []
    for d in DELTAS[:1]:
        v = list(v0)
        ch, st = L.ripple(v, {x: v0[x] + d})
        tou = L.touched_atoms(v, av0, ch)
        at_union |= set(tou)
        # recompute equation values only for equations touching a changed atom
        cand = L.eqs_of_atoms(tou)
        avn = list(av0)
        for a, nv in tou.items():
            avn[a] = nv
        eqs = set(i for i in cand if L.eq_value(i, avn) != base_eq[i])
        eq_union |= eqs
        per_delta.append((len(ch), len(tou), len(eqs)))
    res[x] = {
        'atoms': frozenset(at_union),
        'eqs': frozenset(eq_union),
        'nchanged': per_delta[0][0],
    }
    if x % 5000 == 0:
        print(f'  {x}/{L.NVARS}  {time.time()-t0:.0f}s', file=sys.stderr)

print(f'census done in {time.time()-t0:.0f}s')

# ---- comparison against the syntactic footprint -------------------------------
gain = 0
strict = 0
empty_eff = 0
for x in range(L.NVARS):
    syn = L.var_eqs.get(x, set())
    eff = res[x]['eqs']
    if len(eff) < len(syn):
        gain += 1
        strict += len(syn) - len(eff)
    if not eff:
        empty_eff += 1
print(f'variables whose EFFECTIVE eq-footprint is strictly smaller than the SYNTACTIC one: {gain}')
print(f'total equations removed across all variables: {strict}')
print(f'variables with EMPTY effective footprint (pure no-ops): {empty_eff}')
hist = collections.Counter(len(res[x]['eqs']) for x in range(L.NVARS))
print('effective footprint size histogram (size:count) for sizes <= 30:')
for k in sorted(hist):
    if k <= 30:
        print(f'   {k}: {hist[k]}')

pickle.dump({'res': res, 'base_eq': base_eq, 'fails0': fails0, 'path': base_path},
            open(f'census{which}.pkl', 'wb'))
print(f'wrote census{which}.pkl')
