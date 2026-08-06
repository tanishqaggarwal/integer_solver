"""S10 step 29: verify precisely that x_28730 += k*p closes, with the seed FORBIDDEN,
and report exactly what the resulting state scores and what K2 is."""
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
FORBID = {28730, 4432}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
base_av = L.all_atom_values(base)
BASE_NZ = set(a for a in range(L.NA) if base_av[a])


def bad_after(v, changed):
    t = L.touched_atoms(v, base_av, changed)
    return sorted(set(a for a, val in t.items() if val != 0) - BASE_NZ)


for k in (1, 2, -5):
    d = k * P
    v = list(base)
    ch, _ = L.ripple(v, {28730: base[28730] + d, 4432: base[4432] + d}, block=BLOCK)
    bad = bad_after(v, ch)
    print(f'\nk={k}: initial collateral {bad}')
    frontier = [(v, ch, bad, ())]
    done = None
    for depth in range(6):
        nxt = []
        for vv, cc, bb, path in frontier:
            for a in bb:
                for u in sorted(L.avars[a]):
                    if (a, u) in path or u in FORBID:
                        continue
                    nv = T.solve_lin(a, u, vv)
                    if nv is None or nv == vv[u]:
                        continue
                    w = list(vv)
                    try:
                        ch2, _ = L.ripple(w, {u: nv}, block=BLOCK)
                    except Exception:
                        continue
                    allch = dict(cc); allch.update(ch2)
                    b2 = bad_after(w, allch)
                    if not b2:
                        done = (w, path + ((a, u),)); break
                    nxt.append((w, allch, b2, path + ((a, u),)))
                if done: break
            if done: break
        if done: break
        nxt.sort(key=lambda t: len(t[2]))
        frontier = nxt[:150]
        if not frontier: break
    if done:
        w, path = done
        av = L.all_atom_values(w)
        fail = L.failing_eqs(av)
        print(f'  CLOSED via {path}')
        print(f'  nonzero atoms = {[x for x in range(L.NA) if av[x]]}')
        print(f'  failing={len(fail)} score={L.NEQ-len(fail)}')
        print(f'  K2 before={base[28730] % P}')
        print(f'  K2 after ={w[28730] % P}   (unchanged: {w[28730] % P == base[28730] % P})')
        print(f'  a22230 before={base_av[22230]}')
        print(f'  a22230 after ={av[22230]}')
        print(f'  delta/p = {(av[22230]-base_av[22230])//P if (av[22230]-base_av[22230]) % P == 0 else "not a multiple of p"}')
    else:
        print('  did not close')
