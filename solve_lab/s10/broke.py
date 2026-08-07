import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
def pr(a, n=150):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
print('--- atoms broken by moving x_14623 (a21617 route) ---')
for a in (25676, 33796, 42245):
    print(f'a{a} ({len(L.atom2eq[a])} eqs) {"->x_"+str(atom_out[a][1]) if a in atom_out else "CHECK"}: {pr(a)}')
print('\n--- atoms broken by moving x_14853 (a29539 route) ---')
for a in (19297, 19299, 30984, 36185, 40812):
    print(f'a{a} ({len(L.atom2eq[a])} eqs) {"->x_"+str(atom_out[a][1]) if a in atom_out else "CHECK"}: {pr(a)}')
print('\n--- consumers of the two broadcast-constant inputs ---')
for u in (14623, 14853):
    print(f'x_{u} (val {str(v[u])[:30]}, {v[u].bit_length()} bits):')
    for a in sorted(L.var_atoms[u]):
        print(f'   a{a:<6} ({len(L.atom2eq[a]):>2} eqs) '
              f'{"->x_"+str(atom_out[a][1]) if a in atom_out else "CHECK":>12}: {pr(a, 110)}')
