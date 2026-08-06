"""WR step 4: the exact structure of a37694's 12 equations (+ a39417's 1)."""
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
F = W.F_WIRE
b2 = list(base); F.fwd(b2)
av0 = L.all_atom_values(b2)

E = sorted(set(L.atom2eq[37694]) | set(L.atom2eq[39417]))
print(f'{len(E)} equations carrying a37694 or a39417\n')
allat = collections.Counter()
for e in E:
    m, sq, co = L.eq_atoms[e]
    print(f'eq{e}: mult={m} sq={sq} natoms={len(co)}')
    for a, c in sorted(co.items()):
        tag = 'GATE' if a in F.atom_out else 'chk '
        wv = 'WIRE-ONLY' if L.avars[a] and all(u in WIRE for u in L.avars[a]) else ''
        print(f'    {c:>6} * a{a:<6} {tag} neq={len(L.atom2eq[a]):<3} '
              f'nv={len(L.avars[a]):<3} {wv} {L.atom_src[a][:70]}')
        allat[a] += 1
    print()
print('atom frequency across those equations:')
for a, n in allat.most_common():
    tag = 'GATE' if a in F.atom_out else 'CHECK'
    print(f'   a{a:<6} in {n} of them   {tag}  total_neq={len(L.atom2eq[a])}')
