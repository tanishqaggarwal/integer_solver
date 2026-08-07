"""Full generator scan at the 39,026 placement over ALL variables (blocked gates),
including the output variables of the deliberately broken gates."""
import json, time, sys
import dlib as L
import engine2 as E

w = L.load('../best/new_instance_partial_39026.json')
st0 = E.St(w)
NZ = st0.nz()
Eq = set()
for a in NZ:
    Eq |= set(L.atom2eq.get(a, {}))
atoms = set()
for i in Eq:
    atoms |= set(L.eq_atoms[i][2])
COLS = [a for a in sorted(atoms) if set(L.atom2eq.get(a, {})) <= Eq]
BLOCK = set(a for a in NZ if a in L.atom_out)
st0 = E.St(w, block=BLOCK)
print('score', st0.score, 'block', sorted(BLOCK), 'COLS', COLS)
base = [st0.av[a] for a in COLS]

DELTAS = [1]
gens = []
t0 = time.time()
for u in range(L.NVARS):
    for d in DELTAS:
        st = st0.clone()
        st.apply({u: st.v[u] + d})
        bad = False
        for a in range(L.NA):
            if st.av[a] != st0.av[a] and a not in atoms:
                bad = True
                break
        if bad:
            continue
        dd = [st.av[a] - base[i] for i, a in enumerate(COLS)]
        if any(dd):
            gens.append((u, d, dd))
    if u % 5000 == 0:
        print('  ', u, f'{time.time()-t0:.0f}s gens={len(gens)}', flush=True)
print('cost-free generators over ALL variables:', len(gens))
for u, d, dd in gens:
    print('  x_%-7d d=%+d %s' % (u, d, [(COLS[i], str(x)[:22]) for i, x in enumerate(dd) if x]))
json.dump({'cols': COLS, 'base': [str(x) for x in base],
           'gens': [[u, [str(x) for x in dd]] for u, d, dd in gens]}, open('gens26b.json', 'w'))
