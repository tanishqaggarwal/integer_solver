"""jm step 4: frame-2 repair engine.

After a move, walk the broken CHECK atoms and try to zero each of them through a
FREE input (1-level: the atom's own free variable; 2-level: a variable of the
atom that is a gate output, solved back through its definer to a free input).
Accept any strictly-better (out12, -#broken) state.  Repeat to a fixed point.
"""
import os, sys, time, json, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
P = J.P
w = J.base_state()
R0 = J.resid(w)

FORBID = set()


def state(v):
    c, s, f, av = J.cost(v)
    nz = [a for a in range(L.NA) if av[a] and a not in J.SS]
    return c, s, nz, av


def congr(v):
    """which relaxations this state has achieved, relative to base"""
    av = L.all_atom_values(v)
    return ((v[7068] - v[2099]) % P, (v[28730]) % P, (v[14853] - v[1308]) % P,
            (v[24548] - v[25442]) % P)


def free_targets(a, v):
    """candidate (free input u, new value) that zero atom a."""
    out = []
    for u in sorted(L.avars[a]):
        if u in FORBID:
            continue
        if u in J.FREESET:
            nv = T.solve_lin(a, u, v)
            if nv is not None and nv != v[u]:
                out.append((u, nv))
        else:
            tgt = T.solve_lin(a, u, v)
            if tgt is None or tgt == v[u]:
                continue
            d = J.definer.get(u)
            if d is None:
                continue
            vv = list(v)
            vv[u] = tgt
            for z in sorted(L.avars[d]):
                if z == u or z not in J.FREESET or z in FORBID:
                    continue
                nv = T.solve_lin(d, z, vv)
                if nv is not None and nv != v[z]:
                    out.append((z, nv))
    return out


def repair(v, maxit=25, verbose=True, keep=None):
    """keep: callable(v)->bool, a predicate the repaired state must preserve."""
    c, s, nz, av = state(v)
    best = (-c, -len(nz))
    for it in range(maxit):
        got = None
        for a in nz:
            for u, nv in free_targets(a, v):
                tr = list(v)
                tr[u] = nv
                J.fwd2(tr, 2)
                if keep is not None and not keep(tr):
                    continue
                c2, s2, nz2, av2 = state(tr)
                k = (-c2, -len(nz2))
                if k > best:
                    got = (a, u, tr, c2, s2, nz2, k)
                    break
            if got:
                break
        if not got:
            break
        a, u, v, c, s, nz, best = got
        if verbose:
            print(f'      repair it{it}: a{a} via x_{u} -> out12={c} score={s} '
                  f'broken={nz}', flush=True)
    return v, c, s, nz


if __name__ == '__main__':
    print('=== J2 with correct tracking order ===')
    v = list(w)
    v[28730] += 1000003
    J.fwd2(v, 2)
    v[24548] += v[25442] - w[25442]
    J.fwd2(v, 2)
    c, s, nz, av = state(v)
    print(f'  after 2 trackers: out12={c} score={s} broken={nz}')
    print(f'  x_27522 moved? {v[27522] != w[27522]}  x_14623={v[14623] == w[14623]} '
          f'a21617={av[21617] != 0} a21617 mod p = {av[21617] % P != 0}')
    v2 = list(v)
    v2[14623] += v2[27522] - w[27522]
    J.fwd2(v2, 2)
    c, s, nz, av2 = state(v2)
    print(f'  + track x_27522 with x_14623: out12={c} score={s} broken={nz}')

    print('\n=== repair engine on J2 ===')
    vr, c, s, nz = repair(list(v))
    print(f'  J2 repaired: out12={c} score={s} broken={nz}')
    print(f'  x_28730 mod p changed: {(vr[28730] - w[28730]) % P != 0}')

    print('\n=== repair engine on x_6418 move ===')
    v = list(w); v[6418] += 1000003; J.fwd2(v, 2)
    vr, c, s, nz = repair(list(v))
    print(f'  x_6418 repaired: out12={c} score={s} broken={nz}')
    print(f'  R_A changed: {J.resid(vr)[0] != R0[0]}')

    print('\n=== repair engine on J1 (x_7068 + x_14853) ===')
    v = list(w); v[7068] += 1000003; J.fwd2(v, 2)
    v[14853] += v[1308] - w[1308]; J.fwd2(v, 2)
    vr, c, s, nz = repair(list(v))
    print(f'  J1 repaired: out12={c} score={s} broken={nz}')
    print(f'  R_A changed: {J.resid(vr)[0] != R0[0]}')
