"""S10 step 21: are the handles x_27676 = x_23333*x_6504 and x_7574 = x_8143*x_26658
dormant products that can be ACTIVATED?  If so x_31861 and x_14865 become free knobs
on the MUX branch, giving 4 knobs for the 4 conditions."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(HERE, 'mux_on.json'))

for name, (a, b) in (('x_27676', (23333, 6504)), ('x_7574', (8143, 26658))):
    print(f'\n=== {name} = x_{a} * x_{b} ===')
    for u in (a, b):
        d = L.definer.get(u)
        print(f'  x_{u:<7} val={v[u]!s:<28} free={u not in L.definer} '
              f'natoms={len(L.var_atoms[u]):<3} neqs={len(L.var_eqs[u])}')
        if d is not None:
            print(f'      definer a{d}: {L.atom_src[d][:130]}')
        for at in sorted(L.var_atoms[u]):
            print(f'        in a{at:<6} (neq={len(L.atom2eq.get(at,{})):<3}) '
                  f'{L.atom_src[at][:95]}')

print('\n=== can we set them nonzero?  perturb each factor and count damage ===')
base_av = L.all_atom_values(v)
BASE_NZ = set(x for x in range(L.NA) if base_av[x])
for u in (23333, 6504, 8143, 26658):
    w = list(v)
    L.ripple(w, {u: 1 if v[u] == 0 else v[u] + 1})
    wav = L.all_atom_values(w)
    nz = set(x for x in range(L.NA) if wav[x])
    print(f'  x_{u:<7} -> extra_nz={sorted(nz-BASE_NZ)} gone={sorted(BASE_NZ-nz)} '
          f'x_27676={w[27676]} x_7574={w[7574]}')
