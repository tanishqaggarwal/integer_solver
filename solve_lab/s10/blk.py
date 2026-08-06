import os, sys
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
def pr(a, n=170):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'
w = L.load(os.path.join(HERE, 'jm_azero00_39017.json'))
d = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
FREE = set(ad.FREE)
for a in (688, 1618, 40608, 687, 1617):
    print(f'\na{a} ({len(L.atom2eq[a])} eqs) '
          f'{"->x_"+str(L.atom_out[a][1]) if a in L.atom_out else "CHECK"}')
    print(f'   {pr(a)}')
    for x in sorted(set(L.avars[a])):
        print(f'      x_{x:<7} {"FREE" if x in FREE else "gate a"+str(L.definer.get(x)):<12}'
              f' 39017val={str(w[x])[:20]:<22} delivered={str(d[x])[:20]:<22}'
              f' equal_mod_p={w[x]%P == d[x]%P}')
