"""WR step 6: at w=1, try to zero the 7 residual atoms directly through their
handles (which are unquantised once the wire is 1)."""
import os, sys, collections, itertools
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
P = ad.P

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
WIRE = W.wire_of(base)
F = W.F_WIRE
b2 = list(base); F.fwd(b2)
WV = int(sys.argv[1]) if len(sys.argv) > 1 else 1
v = list(b2)
for u in WIRE:
    v[u] = WV
F.fwd(v, rounds=10)
av, nz, fail, sc = F.report(v, f'w={WV} start')

RES = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
print('\nresidual atoms and the variables available in each:')
for a in RES:
    print(f'  a{a}: {L.atom_src[a][:90]}')
    for u in sorted(L.avars[a]):
        kind = ('WIRE' if u in set(WIRE) else
                ('DETACHED-free' if u in F.detach else
                 ('free' if u not in F.definer else 'gate')))
        print(f'      x_{u:<7} {kind:<14} used_in_atoms={len(L.var_atoms[u])} '
              f'val_bits={abs(v[u]).bit_length()}  solvable={T.solve_lin(a, u, v) is not None}')

# greedy: repeatedly pick the (atom, var) move that most improves the potential
def pot(vv):
    a2 = L.all_atom_values(vv)
    nz2 = [a for a in range(L.NA) if a2[a]]
    f2 = L.failing_eqs(a2)
    return (L.NEQ - len(f2), -len(nz2)), a2, nz2

cur, av, nz = pot(v)
print(f'\ngreedy from {cur}')
for it in range(40):
    best = None
    for a in list(nz):
        for u in sorted(L.avars[a]):
            if u in W.F_WIRE.detach and u == 26064:
                continue
            tgt = T.solve_lin(a, u, v)
            if tgt is None or tgt == v[u]:
                continue
            tr = list(v); tr[u] = tgt
            F.fwd(tr, rounds=6)
            p2, av2, nz2 = pot(tr)
            if best is None or p2 > best[0]:
                best = (p2, a, u, tr, av2, nz2)
    if best is None or best[0] <= cur:
        print(f'  it{it}: stuck at {cur}; nonzero {sorted(nz)}')
        break
    p2, a, u, tr, av2, nz2 = best
    print(f'  it{it}: a{a} via x_{u}: {cur} -> {p2}', flush=True)
    v, cur, av, nz = tr, p2, av2, nz2
print(f'FINAL score {cur[0]}  nonzero {sorted(nz)}')
T.save(v, os.path.join(HERE, f'wr_rep1_{WV}_{cur[0]}.json'))
