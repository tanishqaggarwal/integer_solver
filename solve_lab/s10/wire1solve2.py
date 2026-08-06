"""S10 step 44: wire=1 solve, routing every fix THROUGH the solo handle.

To zero atom `a` we pick a gate-output t occurring in it, compute the value t_new
that makes a vanish, then set the solo free handle h of t's own defining atom so
that the forward evaluation reproduces t_new.  Because h occurs in exactly one
atom, this costs nothing anywhere else -- and with the wire at 1 the handle is
unquantised, so it always has an exact integer solution.
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
ROOTATOM = 37694
BLOCK = {ROOTATOM}

SOLO = collections.defaultdict(list)
for u in range(L.NVARS):
    if u not in L.definer and len(L.var_atoms[u]) == 1:
        SOLO[L.var_atoms[u][0]].append(u)


def fwd_block(v, rounds=4):
    for _ in range(rounds):
        for u in ad.ORDER:
            a = L.definer[u]
            if a in BLOCK:
                continue
            nv = T.solve_lin(a, u, v)
            if nv is not None:
                v[u] = nv
    return v


def status(v):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    return av, nz, L.failing_eqs(av)


def handle_moves(v, a):
    """[(handle, value, via_t)] making atom a vanish through a solo handle."""
    out = []
    for t in sorted(L.avars[a]):
        r = T.lin_parts(a, t, v)
        if r is None:
            continue
        c, rest = r
        if c == 0 or rest % c:
            continue
        t_new = -rest // c
        if t_new == v[t]:
            continue
        d = L.definer.get(t)
        if d is None:                    # t itself is a free input
            out.append((t, t_new, t))
            continue
        w = list(v); w[t] = t_new
        for h in SOLO.get(d, []):
            nv = T.solve_lin(d, h, w)
            if nv is not None and nv != v[h]:
                out.append((h, nv, t))
    return out


v = L.load(os.path.join(HERE, 'wire1_state.json'))
av, nz, fail = status(v)
print(f'start: nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)

for it in range(40):
    av, nz, fail = status(v)
    todo = [a for a in nz if a != ROOTATOM]
    print(f'iter {it}: nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)
    if not todo:
        print('*** only the wire pin a37694 remains ***', flush=True)
        break
    best = None
    for a in todo:
        for h, val, t in handle_moves(v, a):
            w = list(v); w[h] = val
            fwd_block(w)
            av2, nz2, f2 = status(w)
            key = (len(nz2), len(f2))
            if best is None or key < best[0]:
                best = (key, a, h, val, t, w)
    if best is None:
        print('  no handle move available', flush=True); break
    key, a, h, val, t, w = best
    if key[0] >= len(nz):
        print(f'  best move (a{a} via x_{h}) does not reduce: {key}', flush=True); break
    print(f'  a{a}: handle x_{h} (via x_{t}) -> nz={key[0]} failing={key[1]}', flush=True)
    v = w

av, nz, fail = status(v)
print(f'\nFINAL nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}')
T.save(v, os.path.join(HERE, 'wire1_solved2.json'))
if nz == [ROOTATOM] or not nz:
    print('\n=== the wire pin is the only obstruction; its 12 equations: ===')
    for i in sorted(L.atom2eq.get(ROOTATOM, {})):
        m, sq, co = L.eq_atoms[i]
        print(f'  eq {i:<6} coeff_of_a37694={co[ROOTATOM]:<5} n_atoms={len(co)} '
              f'value={L.eq_value(i, av)!s:.30}')
