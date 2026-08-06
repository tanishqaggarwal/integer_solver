"""S10 step 8: actually construct the MUX branch and measure it in equation space."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')

print('=== how x_21279 is switched on ===')
for u in (9062, 20434, 4287, 2081):
    d = L.definer.get(u)
    print(f'  x_{u:<6} def_atom={str(d):<8} '
          f'{L.atom_src[d][:100] if d is not None else "FREE INPUT"}')

v = L.load(BEST)
print('  current: x_9062=%s x_20434=%s x_4287=%s x_2081=%s x_21279=%s'
      % (v[9062], v[20434], v[4287], v[2081], v[21279]))

def measure(v, tag):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'  [{tag}] nonzero atoms={len(nz)} failing={len(fail)} score={L.NEQ-len(fail)}')
    return av, nz, fail

v = L.load(BEST)
measure(v, 'base')

# switch the MUX on
ch, _ = L.ripple(v, {4287: 1})
print(f'\nafter x_4287:=1 ripple: {len(ch)} vars changed; x_21279={v[21279]} x_7075={v[7075]}')
av, nz, fail = measure(v, 'mux on')
for a in nz:
    print(f'    a{a:<6} neq={len(L.atom2eq.get(a,{})):<3} {L.atom_src[a][:95]}')
json.dump({f'x_{i}': v[i] for i in range(L.NVARS) if v[i] != 0},
          open(os.path.join(HERE, 'mux_on.json'), 'w'))
print('saved mux_on.json')
