"""S10 step 32: Gauss-Seidel on the CHECK system in the forward-eval frame.

Forward eval leaves exactly 6 failing checks, every one containing a free input:
  a7930   x_24548 == x_25442 (mod p)      x_24548 FREE
  a29539  x_14853 == x_1308  (mod p)      x_14853 FREE
  a35759  x_9118  == 0       (mod p)      x_9118  FREE
  a35760  x_8731  == 0       (mod p)      x_8731  FREE
  a40826, a41512                          big checks, 1 equation each

Iterate: forward-eval -> solve each failing check for one of its free inputs ->
forward-eval again.  If this converges, every atom is zero and the instance is solved.
"""
import os, sys, json, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
topo = list(L.topo)
cyc = [x for x in definer if x not in set(topo)]
FREE = set(range(L.NVARS)) - set(definer)


def fwd(v, rounds=3):
    for _ in range(rounds):
        for seq in (topo, cyc):
            for u in seq:
                a = definer[u]
                nv = T.solve_lin(a, u, v)
                if nv is not None:
                    v[u] = nv
    return v


def status(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    return av, nz, fail


v = L.load(os.path.join(HERE, 'forward_state.json'))
av, nz, fail = status(v)
print(f'start: nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)

best = (L.NEQ - len(fail), list(v))
history = []
for it in range(40):
    av, nz, fail = status(v)
    sc = L.NEQ - len(fail)
    if sc > best[0]:
        best = (sc, list(v))
    print(f'\niter {it}: nz={nz} failing={len(fail)} score={sc}', flush=True)
    if not nz:
        print('*** ALL ATOMS ZERO -- FULL SOLUTION ***', flush=True)
        break
    key = tuple(nz)
    if history.count(key) >= 2:
        print('  cycling; stop', flush=True)
        break
    history.append(key)
    moved = False
    for a in nz:
        # prefer a FREE input; among those prefer the one in fewest atoms
        cands = sorted((u for u in L.avars[a] if u in FREE),
                       key=lambda u: (len(L.var_atoms[u]), len(L.var_eqs[u])))
        for u in cands:
            nv = T.solve_lin(a, u, v)
            if nv is None or nv == v[u]:
                continue
            v[u] = nv
            fwd(v)
            av2 = L.all_atom_values(v)
            print(f'   a{a}: set x_{u} (free, {len(L.var_atoms[u])} atoms) -> '
                  f'nz now {[x for x in range(L.NA) if av2[x]]}', flush=True)
            moved = True
            break
        if moved:
            break
    if not moved:
        print('  no free-input move available', flush=True)
        break

av, nz, fail = status(v)
print(f'\nFINAL: nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}')
print(f'BEST seen: {best[0]}')
T.save(best[1], os.path.join(HERE, 'gs_best.json'))
