"""S10 step 19: attack the two congruences with a STRONGER repair than the ripple.

The ripple only repairs an atom through its canonical output var with exact
division.  repairD.py showed that explicitly choosing a *handle* variable fixes
atoms the ripple gives up on.  So re-run the two forbidden moves with a
best-first repair over ALL variables of every broken atom.

  move A: x_28730 += 1        (kills congruence 2 if 7930/41512 can be closed)
  move B: x_7068  += 1        (kills congruence 1 if 29539/40826 can be closed)
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BLOCK = set(NZ) | {22231}
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
base = L.load(BEST)
BASE_NZ = set(a for a in range(L.NA) if L.all_atom_values(base)[a])


def repair(v, rounds=8, verbose=True):
    """Best-first: repeatedly close the broken atom/var pair that minimises the
    number of nonzero atoms outside the residual set."""
    seen = set()
    for r in range(rounds):
        av = L.all_atom_values(v)
        bad = [a for a in range(L.NA) if av[a] and a not in BASE_NZ]
        if not bad:
            return v, []
        best = None
        for a in bad:
            for u in sorted(L.avars[a]):
                if (a, u) in seen:
                    continue
                nv = T.solve_lin(a, u, v)
                if nv is None or nv == v[u]:
                    continue
                w = list(v)
                try:
                    L.ripple(w, {u: nv}, block=BLOCK)
                except Exception:
                    continue
                wav = L.all_atom_values(w)
                nbad = len([x for x in range(L.NA) if wav[x] and x not in BASE_NZ])
                key = (nbad, len(L.var_eqs[u]))
                if best is None or key < best[0]:
                    best = (key, a, u, nv, w)
        if best is None:
            break
        key, a, u, nv, w = best
        seen.add((a, u))
        if verbose:
            print(f'    round {r}: close a{a} via x_{u} -> outstanding={key[0]}')
        v = w
        if key[0] == 0:
            return v, []
    av = L.all_atom_values(v)
    return v, [a for a in range(L.NA) if av[a] and a not in BASE_NZ]


for tag, seeds in (('A  x_28730 += 1', {28730: base[28730] + 1, 4432: base[4432] + 1}),
                   ('B  x_7068  += 1', {7068: base[7068] + 1})):
    print(f'\n===== {tag} =====')
    v = list(base)
    L.ripple(v, seeds, block=BLOCK)
    av = L.all_atom_values(v)
    print('  immediately broken:', [a for a in range(L.NA) if av[a] and a not in BASE_NZ])
    v, left = repair(v)
    av = L.all_atom_values(v)
    fail = L.failing_eqs(av)
    print(f'  after repair: outstanding={left} failing={len(fail)} score={L.NEQ-len(fail)}')
    print(f'    D%p  = {(v[7068]-v[2099]) % P}')
    print(f'    K2   = {v[28730] % P}')
    if not left:
        T.save(v, os.path.join(HERE, f'congr_{tag.split()[0]}.json'))
        print(f'    *** congruence moved with NO collateral -> saved')
