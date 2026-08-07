import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
J = json.load(open(os.path.join(HERE, 'cone.json')))
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
cone = set(J['cone'])
v = L.load(os.path.join(HERE, 'forward_state.json'))
av = L.all_atom_values(v)

def term(m, c):
    if not m: return f'{c}'
    s = '*'.join(f'x_{w}' for w in m)
    if c == 1: return s
    if c == -1: return '-' + s
    return f'{c}*{s}'

def pretty(a):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    return ' + '.join(term(m, c) for m, c in ts).replace('+ -', '- ')

print(f'{"var":>9} {"definer":>9} {"value (trunc)":>30}')
for t in sorted(cone):
    d = definer.get(t)
    print(f'x_{t:<7} {"a"+str(d) if d is not None else "FREE":>9} '
          f'{str(v[t])[:30]:>30}  bits={v[t].bit_length() if isinstance(v[t],int) else "?"}')

print('\n=== the 24 atoms of the cone ===')
inside = [a for a in range(L.NA) if L.avars[a] and set(L.avars[a]) <= cone]
for a in sorted(inside):
    tag = f'-> x_{atom_out[a]}' if a in atom_out else '[CHECK]'
    print(f'a{a:<6} {tag:>10} eqs={len(L.atom2eq[a]):>3} val={str(av[a])[:24]:>26}')
    print(f'          {pretty(a)[:300]}')
