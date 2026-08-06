"""S10 step 90: the exact variable footprint of the seven residual atoms."""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]

def term(m, c):
    if not m: return f'{c}'
    s = '*'.join(f'x_{w}' for w in m)
    return s if c == 1 else ('-' + s if c == -1 else f'{c}*{s}')
def pretty(a):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    return ' + '.join(term(m, c) for m, c in ts).replace('+ -', '- ')

for a in SEVEN:
    print(f'\n===== a{a}  ({len(L.atom2eq[a])} eqs)  defines '
          f'{"x_"+str(atom_out[a]) if a in atom_out else "(check)"} =====')
    print('   ', pretty(a))
    for w in sorted(set(L.avars[a])):
        others = sorted(set(L.var_atoms[w]) - {a})
        eqs = sum(len(L.atom2eq[o]) for o in others)
        tag = 'SOLO HANDLE' if not others else f'{len(others)} other atoms, {eqs} eqs'
        print(f'      x_{w:<7} = {str(v[w])[:26]:<28} bits={v[w].bit_length():<5} {tag}')
        for o in others:
            print(f'            a{o:<6} eqs={len(L.atom2eq[o]):>3} '
                  f'{"defines x_"+str(atom_out[o]) if o in atom_out else "CHECK"}  '
                  f'{"ZERO" if av[o]==0 else "NONZERO"}  vars={len(set(L.avars[o]))}')
