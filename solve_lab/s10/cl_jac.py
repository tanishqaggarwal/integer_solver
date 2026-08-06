"""CL step 5: the exact mod-p Jacobian of every CHECK w.r.t. the free inputs that can
move the two failing checks, and the left-kernel obstruction certificate."""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, fwdad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
vm = [x % P for x in v0]
av0 = L.all_atom_values(v0)
CHECKS = [a for a in range(L.NA) if a not in atom_out]
BAD = [a for a in CHECKS if av0[a]]
print(f'{len(CHECKS)} checks, failing {BAD}')

U = set()
for a in BAD:
    U |= set(ad.grad(a, vm))
U = sorted(U - {2081, 4287})
print(f'free inputs able to move a failing check: {len(U)}')

t0 = time.time()
cols = {}
for i, u in enumerate(U):
    cols[u] = fwdad.jac_column(u, v0, vm, CHECKS)
    if i % 20 == 0:
        print(f'  col {i}/{len(U)} ({time.time()-t0:.0f}s) support {len(cols[u])}', flush=True)
rows = sorted(set().union(*[set(c) for c in cols.values()]))
print(f'\nJacobian {len(rows)} x {len(U)}   ({time.time()-t0:.0f}s)')
print(f'affected checks currently failing: {[r for r in rows if av0[r]]}')
json.dump({'U': U, 'rows': rows,
           'cols': {str(u): {str(c): str(d) for c, d in cols[u].items()} for u in U},
           'resid': {str(r): str(av0[r] % P) for r in rows}},
          open(os.path.join(HERE, 'cl_jac.json'), 'w'))
print('saved cl_jac.json')
