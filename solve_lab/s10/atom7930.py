"""S10 step 24: atom 7930 is the ONLY thing making congruence (2) binding.

  a7930 = 9367949*(x_24548 - x_25442) - x_7927

If x_28730 (hence x_4432) moves by d, something in {x_24548, x_25442, x_7927}
must absorb it.  Enumerate the structure and test every absorber explicitly.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ) | {22231}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
BASE_NZ = set(a for a in range(L.NA) if L.all_atom_values(base)[a])

for a in (7930, 41512):
    print(f'\na{a}: {L.atom_src[a][:260]}')
    for u in sorted(L.avars[a])[:20]:
        d = L.definer.get(u)
        print(f'   x_{u:<7} val={str(base[u])[:20]:<22} free={u not in L.definer} '
              f'natoms={len(L.var_atoms[u]):<3} neqs={len(L.var_eqs[u])}')
        if d is not None:
            print(f'       def a{d}: {L.atom_src[d][:110]}')

print('\n=== move x_28730 by d, then try every absorber for 7930 / 41512 ===')
for d in (1, P, 7376877):
    v = list(base)
    L.ripple(v, {28730: base[28730] + d, 4432: base[4432] + d}, block=BLOCK)
    av = L.all_atom_values(v)
    broken = sorted(set(x for x in range(L.NA) if av[x]) - BASE_NZ)
    print(f'\n d={"p" if d==P else d}: broken={broken}')
    for a in broken:
        opts = []
        for u in sorted(L.avars[a]):
            nv = T.solve_lin(a, u, v)
            if nv is None or nv == v[u]:
                continue
            w = list(v)
            try: L.ripple(w, {u: nv}, block=BLOCK)
            except Exception: continue
            wav = L.all_atom_values(w)
            extra = sorted(set(x for x in range(L.NA) if wav[x]) - BASE_NZ)
            opts.append((len(extra), u, extra))
        opts.sort()
        for n, u, extra in opts[:4]:
            print(f'    a{a} via x_{u:<7} -> remaining collateral {extra}')
        if not opts:
            print(f'    a{a}: NO absorber')
