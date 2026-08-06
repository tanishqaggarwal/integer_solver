import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from repair2 import score
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
av = L.all_atom_values(v)
print(f'state score {score(v)}; nonzero {[a for a in range(L.NA) if av[a]]}')
for a in [21617, 29539]:
    print(f'\n=== a{a} ({len(L.atom2eq[a])} eqs) val={str(av[a])[:34]}')
    for w in sorted(set(L.avars[a])):
        others = sorted(set(L.var_atoms[w]) - {a})
        print(f'  x_{w:<7} {"FREE" if w in FREE else "def a"+str(definer[w]):<10} '
              f'val={str(v[w])[:22]:<24} bits={v[w].bit_length():<5} '
              f'others={[(o, len(L.atom2eq[o])) for o in others]}')
        for o in others:
            if o in atom_out and len(set(L.avars[o])) <= 4:
                ts = sorted(L.polys[o].items(), key=lambda kv: (len(kv[0]), kv[0]))
                pr = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                                 ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                                  f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                                for m, c in ts).replace('+ -', '- ')
                print(f'        a{o}: {pr[:120]}')
    # divisibility structure
    print('  divisibility:')
    for w in sorted(set(L.avars[a])):
        for c in {abs(cc) for cc in L.polys[a].values()}:
            if c > 1000 and v[w] and v[w] % c == 0:
                print(f'     {c} divides x_{w}')
