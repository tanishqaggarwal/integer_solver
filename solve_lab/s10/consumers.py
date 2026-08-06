"""S10 step 89: who CONSUMES the four pinned variables?"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)

def term(m, c):
    if not m: return f'{c}'
    s = '*'.join(f'x_{w}' for w in m)
    return s if c == 1 else ('-' + s if c == -1 else f'{c}*{s}')
def pretty(a):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    out = ' + '.join(term(m, c) for m, c in ts).replace('+ -', '- ')
    return out if len(out) < 220 else out[:220] + ' ...'

for t in [7068, 2099, 28730, 642, 9413, 17325, 6418, 9118, 31861]:
    ats = sorted(L.var_atoms[t])
    print(f'\n=== x_{t}: appears in {len(ats)} atoms ===')
    for a in ats:
        role = f'defines x_{atom_out[a]}' if a in atom_out else 'CHECK'
        print(f'  a{a:<6} [{role:>16}] eqs={len(L.atom2eq[a]):>3} '
              f'val={"0" if av[a]==0 else "NONZERO"}')
        print(f'        {pretty(a)}')
