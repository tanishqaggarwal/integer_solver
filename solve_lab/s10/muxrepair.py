"""S10 step 23: repair the collateral atoms on the zeroed-residual MUX state."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ) | {22231}
v = L.load(os.path.join(HERE, 'muxzero.json'))
av = L.all_atom_values(v)
print('start nz:', [a for a in range(L.NA) if av[a]],
      'failing', len(L.failing_eqs(av)))

seen = set()
for rnd in range(20):
    av = L.all_atom_values(v)
    bad = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'round {rnd}: nz={bad} failing={len(fail)} score={L.NEQ-len(fail)}')
    if not bad:
        print('ALL ATOMS ZERO -- FULL SOLVE'); break
    best = None
    for a in bad:
        for u in sorted(L.avars[a]):
            if (a, u) in seen: continue
            nv = T.solve_lin(a, u, v)
            if nv is None or nv == v[u]: continue
            w = list(v)
            try: L.ripple(w, {u: nv}, block=BLOCK)
            except Exception: continue
            wav = L.all_atom_values(w)
            nb = len([x for x in range(L.NA) if wav[x]])
            nf = len(L.failing_eqs(wav))
            key = (nb, nf)
            if best is None or key < best[0]:
                best = (key, a, u, nv, w)
    if best is None:
        print('  no move'); break
    key, a, u, nv, w = best
    seen.add((a, u))
    if key[0] >= len(bad) and rnd > 3:
        print(f'  best move a{a} via x_{u} gives nz={key[0]} fail={key[1]}; no progress')
        break
    print(f'  -> close a{a} via x_{u}: nz={key[0]} failing={key[1]}')
    v = w
av = L.all_atom_values(v)
fail = L.failing_eqs(av)
print(f'\nFINAL: nz={[a for a in range(L.NA) if av[a]]} failing={len(fail)} '
      f'score={L.NEQ-len(fail)}')
T.save(v, os.path.join(HERE, 'muxrepair.json'))
