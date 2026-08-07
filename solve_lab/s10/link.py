"""S11 step 55: is the a7930 link  x_25442 == x_12553 (mod p)  frame-independent?

If it is, then with the selectors on the chain forces C3 == K1 (mod p), which is
FALSE -- and the instance would have no full solution in that branch.  If the link
is frame-dependent, the infeasibility argument does not close.
"""
import os, sys
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
def pr(a, n=160):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'
print(f'x_25442 is defined by a{L.definer.get(25442)}')
d = L.definer.get(25442)
if d is not None: print(f'   a{d}: {pr(d)}')
print(f'x_25442 appears in atoms {sorted(L.var_atoms[25442])}')
for a in sorted(L.var_atoms[25442]): print(f'   a{a} ({len(L.atom2eq[a])} eqs): {pr(a)}')
print(f'\nx_13458 (the a3578 slack) is defined by a{L.definer.get(13458)}')
d2 = L.definer.get(13458)
if d2 is not None: print(f'   a{d2}: {pr(d2)}')

for name, path in [('delivered', os.path.join(LAB,'best','new_instance_partial_39026.json')),
                   ('canonical 39,009', os.path.join(HERE,'mod9118_0.json'))]:
    v = L.load(path)
    av = L.all_atom_values(v)
    print(f'\n=== {name} ===')
    print(f'  a7930 = {"ZERO" if av[7930]==0 else "nonzero"}   '
          f'a3578 = {"ZERO" if av[3578]==0 else "nonzero"}   '
          f'a21617 = {"ZERO" if av[21617]==0 else "nonzero"}   '
          f'a31672 = {"ZERO" if av[31672]==0 else "nonzero"}')
    print(f'  x_25442 mod p == x_12553 mod p ? {v[25442]%P == v[12553]%P}')
    print(f'  x_25442 mod p == x_24548 mod p ? {v[25442]%P == v[24548]%P}')
    print(f'  x_2081 = {v[2081]}, x_24601 = {v[24601]}')
    print(f'  x_25442 mod p = {str(v[25442]%P)[:44]}')
    print(f'  x_12553 mod p = {str(v[12553]%P)[:44]}')
