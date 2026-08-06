"""WR step 3: exact accounting of the w!=p branch."""
import os, sys, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
P = ad.P

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
WIRE = W.wire_of(base)
F = W.F_WIRE
b2 = list(base); F.fwd(b2)
av0, nz0, fail0, sc0 = F.report(b2, 'baseline')

v = list(b2)
for u in WIRE:
    v[u] = 1
F.fwd(v, rounds=10)
av1, nz1, fail1, sc1 = F.report(v, 'w=1')

print()
for a in sorted(nz1):
    print(f'  a{a:<6} neq={len(L.atom2eq[a]):<3} val_bits={abs(av1[a]).bit_length():<5} '
          f'{L.atom_src[a][:110]}')

E7694 = set(L.atom2eq[37694]); E9417 = set(L.atom2eq[39417])
BASE7 = set(fail0)
F1 = set(fail1)
print(f'\neqs(a37694)={len(E7694)}  eqs(a39417)={len(E9417)}  '
      f'overlap={len(E7694 & E9417)}')
print(f'baseline failing (7): {sorted(BASE7)}')
print(f'w=1 failing ({len(F1)}): {sorted(F1)}')
print(f'  new failures: {sorted(F1 - BASE7)}   fixed: {sorted(BASE7 - F1)}')
print(f'  new failures inside eqs(37694)+eqs(39417): '
      f'{sorted((F1 - BASE7) & (E7694 | E9417))}')
print(f'  new failures OUTSIDE those: {sorted((F1 - BASE7) - (E7694 | E9417))}')
print(f'  eqs(37694) satisfied anyway: {sorted(E7694 - F1)}')
print(f'  eqs(39417) satisfied anyway: {sorted(E9417 - F1)}')

# which atoms drive the "outside" new failures?
OUT = sorted((F1 - BASE7) - (E7694 | E9417))
for e in OUT[:30]:
    m, sq, co = L.eq_atoms[e]
    nzs = {a: (co[a], av1[a]) for a in co if av1[a]}
    print(f'   eq{e}: sq={sq} atoms_nonzero={ {a: c for a,(c,_) in nzs.items()} }')
