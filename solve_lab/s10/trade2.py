"""S10 step 47: apply the trade using the measured chain coefficients, then look at
exactly which equations remain and what can compensate them."""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

ROOTATOM = 37694
BLOCK = {ROOTATOM}


def fwd_block(v, rounds=4):
    for _ in range(rounds):
        for u in ad.ORDER:
            d = L.definer[u]
            if d in BLOCK:
                continue
            nv = T.solve_lin(d, u, v)
            if nv is not None:
                v[u] = nv
    return v


def rep(v, tag):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'[{tag}] nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)
    return av, nz, fail


v = L.load(os.path.join(HERE, 'wire1_solved2.json'))
av, nz, fail = rep(v, 'start')

# chain coefficients measured in joint.py: d(a7930)/d(x_11052) = -1,
# d(a29539)/d(x_30163) = -1
for atom, knob in ((7930, 11052), (29539, 30163)):
    av = L.all_atom_values(v)
    w = list(v); w[knob] += 1; fwd_block(w)
    d = L.all_atom_values(w)[atom] - av[atom]
    print(f'  d(a{atom})/d(x_{knob}) = {d}')
    if d and av[atom] % d == 0:
        v[knob] += -av[atom] // d
        fwd_block(v)
        av, nz, fail = rep(v, f'closed a{atom} via x_{knob}')

T.save(v, os.path.join(HERE, 'trade_out.json'))

av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'\nremaining atoms {nz}; failing equations {len(fail)}: {fail}')

# which atoms in the failing equations are adjustable (free granularity-1 handle)?
SOLO = collections.defaultdict(list)
for u in range(L.NVARS):
    if u not in L.definer and len(L.var_atoms[u]) == 1:
        SOLO[L.var_atoms[u][0]].append(u)
import math
ADJ = {}
for a, us in SOLO.items():
    g = 0
    for u in us:
        r = T.lin_parts(a, u, v)
        if r:
            g = math.gcd(g, abs(r[0]))
    if g == 1:
        ADJ[a] = us
print(f'\nadjustable atoms (granularity 1): {len(ADJ)}')

print('\nper failing equation: adjustable atoms available to compensate')
for i in fail:
    m, sq, co = L.eq_atoms[i]
    hits = [a for a in co if a in ADJ]
    print(f'  eq {i:<6} sq={int(sq)} n_atoms={len(co):<3} adjustable={hits}')
