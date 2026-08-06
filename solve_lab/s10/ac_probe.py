"""S12 step 1: orientation probe -- baseline states, failing atoms, cluster."""
import os, sys, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
print(f'NVARS={L.NVARS} NA={L.NA} NEQ={L.NEQ} nfree={len(FREE)}')
for name in ('mod9118_0.json', os.path.join(LAB,'best','new_instance_partial_39026.json')):
    p = name if os.path.isabs(name) else os.path.join(HERE, name)
    if not os.path.exists(p): print('missing', p); continue
    v = L.load(p)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    checks = [a for a in nz if a not in atom_out]
    fail = L.failing_eqs(av)
    eqs = sorted(L.eqs_of_atoms(nz))
    print(f'\n{os.path.basename(p)}: score {L.NEQ-len(fail)}  nonzero atoms {len(nz)} {nz}')
    print(f'   check atoms nonzero: {checks}')
    print(f'   failing eqs: {fail}')
    print(f'   eqs touched by nonzero atoms: {len(eqs)} {eqs}')
    nzv = sum(1 for i in range(L.NVARS) if v[i])
    nzf = sum(1 for u in FREE if v[u])
    print(f'   nonzero vars {nzv}/{L.NVARS}; nonzero free {nzf}/{len(FREE)}')
