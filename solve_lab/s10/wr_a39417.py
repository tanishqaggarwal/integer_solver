import os, sys, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
P = ad.P
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
WIRE = set(W.wire_of(base))
print('a39417 full src:', L.atom_src[39417])
print('vars:', sorted(L.avars[39417]))
print('in wire?', {u: (u in WIRE) for u in sorted(L.avars[39417])})
print('poly:', L.polys[39417])
print('eqs:', L.atom2eq[39417])
e = list(L.atom2eq[39417])[0]
m, sq, co = L.eq_atoms[e]
print(f'\neq{e}: mult={m} sq={sq} atoms={len(co)}')
F = W.F_WIRE
for a, c in sorted(co.items()):
    tag = 'GATE' if a in F.atom_out else 'chk '
    wo = 'WIRE-ONLY' if L.avars[a] and all(u in WIRE for u in L.avars[a]) else ''
    print(f'  {c:>6} * a{a:<6} {tag} neq={len(L.atom2eq[a]):<3} {wo} {L.atom_src[a][:80]}')

# x_26064 usage
print('\nx_26064 appears in atoms:')
for a in L.var_atoms[26064]:
    tag = 'GATE' if a in F.atom_out else 'chk '
    print(f'  a{a:<6} {tag} neq={len(L.atom2eq[a]):<3} {L.atom_src[a][:90]}')
