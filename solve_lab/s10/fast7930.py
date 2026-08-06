"""S10 step 25 (fast): can atom 7930 be closed while x_28730 moves?

Congruence (2) binds ONLY because moving x_28730 drags x_4432, which breaks
atoms 7930 and 41512.  If both can be closed, the score is 39,027 at once.

Speed: never re-evaluate all 42,267 atoms -- only the atoms touched by the
variables that actually changed (lib.touched_atoms).
"""
import os, sys, json, time
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
base_av = L.all_atom_values(base)
BASE_NZ = set(a for a in range(L.NA) if base_av[a])
print('base nonzero atoms:', sorted(BASE_NZ), flush=True)


def bad_after(v, changed):
    """Atoms nonzero outside the residual set, computed incrementally."""
    touched = L.touched_atoms(v, base_av, changed)
    out = set(a for a, val in touched.items() if val != 0) - BASE_NZ
    # atoms that were nonzero in base and are still nonzero are fine (they're in BASE_NZ)
    return sorted(out)


def apply(v, seeds):
    ch, _ = L.ripple(v, dict(seeds), block=BLOCK)
    return ch


print('\n=== structure ===')
for a in (7930, 41512):
    print(f'a{a}: {L.atom_src[a][:200]}')
    for u in sorted(L.avars[a]):
        print(f'    x_{u:<7} free={u not in L.definer} natoms={len(L.var_atoms[u]):<3} '
              f'neqs={len(L.var_eqs[u]):<4} val={str(base[u])[:22]}')

print('\n=== depth-2 search: move x_28730 by d, then close 7930 and 41512 ===', flush=True)
t0 = time.time()
for d in (1, P):
    v = list(base)
    ch = apply(v, {28730: base[28730] + d, 4432: base[4432] + d})
    bad = bad_after(v, ch)
    print(f'\nd={"p" if d == P else d}: collateral {bad}  ({time.time()-t0:.0f}s)', flush=True)
    frontier = [(v, ch, bad)]
    for depth in range(3):
        nxt = []
        for vv, cc, bb in frontier:
            if not bb:
                continue
            a = bb[0]
            for u in sorted(L.avars[a]):
                nv = T.solve_lin(a, u, vv)
                if nv is None or nv == vv[u]:
                    continue
                w = list(vv)
                try:
                    ch2 = apply(w, {u: nv})
                except Exception:
                    continue
                allch = dict(cc); allch.update(ch2)
                b2 = bad_after(w, allch)
                nxt.append((w, allch, b2))
                if not b2:
                    av = L.all_atom_values(w)
                    fail = L.failing_eqs(av)
                    print(f'  *** depth {depth}: closed via x_{u}; nonzero atoms '
                          f'{[x for x in range(L.NA) if av[x]]} failing={len(fail)} '
                          f'score={L.NEQ-len(fail)}', flush=True)
                    T.save(w, os.path.join(HERE, 'freed_28730.json'))
        nxt.sort(key=lambda t: len(t[2]))
        frontier = nxt[:6]
        print(f'  depth {depth}: {len(nxt)} branches, best collateral '
              f'{frontier[0][2] if frontier else "-"}  ({time.time()-t0:.0f}s)', flush=True)
        if not frontier:
            break
print('done', time.time() - t0)
