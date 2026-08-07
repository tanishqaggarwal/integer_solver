"""S11 step 25: the compensator family -- what are they, and what does each cost?"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))
ESET = set(E)
def pr(a, n=150):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'
for b in (22231, 22232, 22233, 22234, 22235, 19087, 35757, 35756, 35755):
    out = sorted(set(L.atom2eq[b]) - ESET)
    ov = atom_out.get(b)
    print(f'\na{b}  ({len(L.atom2eq[b])} eqs, {len(out)} outside)  '
          f'{"defines x_"+str(ov[1]) if ov else "CHECK"}   val '
          f'{"0" if av[b]==0 else "NONZERO"}')
    print(f'   {pr(b)}')
    if ov:
        t = ov[1]
        others = sorted(set(L.var_atoms[t]) - {b})
        print(f'   x_{t} also appears in: {[(o, len(L.atom2eq[o])) for o in others]}')
        for o in others:
            print(f'      a{o}: {pr(o, 110)}')
