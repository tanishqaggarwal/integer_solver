"""WR step 11: diagnose the 23 failures of the kernel deformation."""
import os, sys, json, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
import wr_rows as R
P = ad.P
WIRE, rows, RE = R.WIRE, R.rows, R.RE
DET = dict(W.F3)
for u in WIRE:
    if u in L.definer:
        DET[u] = L.definer[u]
FW = W.Frame(DET)
WS = set(WIRE)

v = L.load(os.path.join(HERE, 'wr_kern0.json'))
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'score {L.NEQ-len(fail)}  failing {len(fail)}  nonzero atoms {len(nz)}')
print(f'failing eqs: {sorted(fail)}')

# classify nonzero atoms
cls = collections.Counter()
for a in nz:
    vs = L.avars[a]
    if vs and all(u in WS for u in vs):
        cls['wire-identity/copy'] += 1
    elif a in FW.atom_out:
        cls['broken gate'] += 1
    else:
        cls['check'] += 1
print(cls)

print('\nper failing equation, its nonzero atoms:')
drv = collections.Counter()
for e in sorted(fail):
    m, sq, co = L.eq_atoms[e]
    nzs = {a: co[a] for a in co if av[a]}
    kinds = []
    for a in nzs:
        vs = L.avars[a]
        k = ('W' if vs and all(u in WS for u in vs)
             else ('G' if a in FW.atom_out else 'C'))
        kinds.append(f'a{a}{k}')
        drv[a] += 1
    print(f'  eq{e:<6} sq={sq} n_nonzero={len(nzs)} {kinds[:10]}')
print('\natoms driving the failures:')
for a, n in drv.most_common(30):
    vs = L.avars[a]
    k = ('WIRE' if vs and all(u in WS for u in vs)
         else ('GATE' if a in FW.atom_out else 'CHECK'))
    print(f'  a{a:<6} in {n} failing eqs, total {len(L.atom2eq[a])}  {k}  '
          f'{L.atom_src[a][:80]}')
