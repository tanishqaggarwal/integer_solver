"""S11 step 51: realise the 69-atom vector in the FULLY DETACHED frame.

Each colliding pair shares a variable that is DEFINED by an atom already in the
support -- so detaching it costs nothing and gives each atom of the pair its own
free parameter.  Build that frame, see whether the values become independently
settable, and measure what the detachments break OUTSIDE the support.
"""
import os, sys, json, collections, math
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
J = json.load(open(os.path.join(HERE, 'kervec.json')))
SUPP = sorted(J['support'])
SUPPSET = set(SUPP)
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SUPP]))
print(f'support {len(SUPP)} atoms, {len(E)} equations')
gates = [a for a in SUPP if a in L.atom_out]
checks = [a for a in SUPP if a not in L.atom_out]
print(f'  gate atoms {len(gates)}, check atoms {len(checks)}')
DET = {L.atom_out[a][1]: a for a in gates}
print(f'  variables to detach: {len(DET)}')

definer = {t: a for t, a in L.definer.items() if t not in DET}
ORDER = [t for t in ad.ORDER if t not in DET]
FREE = set(t for t in range(L.NVARS) if t not in definer)
def fwd(v, rounds=8):
    for _ in range(rounds):
        for u in ORDER:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None: v[u] = nv
    return v
def score(v): return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
w = list(base); fwd(w)
print(f'\ndelivered witness in this frame: score {score(w)} '
      f'(on-manifold: {score(w) == 39026})')
aw = L.all_atom_values(w)
nz = [a for a in range(L.NA) if aw[a]]
print(f'  nonzero atoms: {nz}')
outside = [a for a in nz if a not in SUPPSET]
print(f'  nonzero OUTSIDE the support: {outside}')

# which atoms outside the support do the detached variables touch?
touched = set()
for t in DET:
    for a in L.var_atoms[t]:
        if a not in SUPPSET: touched.add(a)
eqs_out = set()
for a in touched: eqs_out |= set(L.atom2eq[a])
print(f'\ndetached variables also appear in {len(touched)} atoms OUTSIDE the support')
print(f'  those atoms live in {len(eqs_out)} equations')
print(f'  of which already inside the support equations: {len(eqs_out & set(E))}')
print(f'  NET equations at risk: {len(eqs_out - set(E))}')
sample = sorted(touched)[:20]
print(f'  sample outside atoms: {[(a, len(L.atom2eq[a])) for a in sample]}')
