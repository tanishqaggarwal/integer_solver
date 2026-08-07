"""jm step 5: stronger frame-2 repair engine with a KEEP predicate.

Move generation, in order of preference:
  L1  atom's own free variable            (solve_lin directly)
  L2  atom var that is a gate output -> back-solve its definer to a free input
  L3  same, one level deeper
Every candidate is applied, fwd2'd, and the state measured EXACTLY.  A candidate
is accepted only if it improves (-out12, -#broken atoms) AND keeps the predicate.
"""
import os, sys, time, json, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
P = J.P
W = J.base_state()
R0 = J.resid(W)


def state(v):
    c, s, f, av = J.cost(v)
    nz = [a for a in range(L.NA) if av[a] and a not in J.SS]
    return c, s, nz, av


def _back(u, tgt, v, depth, seen, out, forbid):
    """ways to make x_u == tgt by moving one free input, up to `depth` levels."""
    if u in J.FREESET:
        if u not in forbid and tgt != v[u]:
            out.append((u, tgt))
        return
    if depth <= 0 or u in seen:
        return
    seen = seen | {u}
    d = J.definer.get(u)
    if d is None:
        return
    vv = list(v)
    vv[u] = tgt
    for z in sorted(L.avars[d]):
        if z == u or z in forbid:
            continue
        nv = T.solve_lin(d, z, vv)
        if nv is None or nv == v[z]:
            continue
        _back(z, nv, v, depth - 1, seen, out, forbid)


def gen(a, v, depth=3, forbid=frozenset()):
    out = []
    for u in sorted(L.avars[a]):
        if u in forbid:
            continue
        tgt = T.solve_lin(a, u, v)
        if tgt is None or tgt == v[u]:
            continue
        _back(u, tgt, v, depth, frozenset(), out, forbid)
    # dedupe
    seen = set()
    res = []
    for u, nv in out:
        if (u, nv) in seen:
            continue
        seen.add((u, nv))
        res.append((u, nv))
    return res


def repair(v, keep=None, maxit=40, depth=3, forbid=frozenset(), verbose=True,
           tag=''):
    c, s, nz, av = state(v)
    best = (-c, -len(nz))
    for it in range(maxit):
        got = None
        for a in nz:
            for u, nv in gen(a, v, depth, forbid):
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
            print(f'    {tag} it{it}: a{a} via x_{u} -> out12={c} score={s} '
                  f'broken={nz}', flush=True)
    return v, c, s, nz


# ---- predicates ---------------------------------------------------------
def keep_c1(v):
    """congruence 1 relaxed: C0 = x_7068 - x_2099 moved mod p"""
    return (v[7068] - v[2099]) % P != R0[0]


def keep_c2(v):
    """congruence 2 relaxed: A1 = x_28730 - p*x_9413 moved mod p"""
    return (v[28730] - W[28730]) % P != 0


def keep_both(v):
    return keep_c1(v) and keep_c2(v)


if __name__ == '__main__':
    t0 = time.time()
    print('=== congruence 2, cheapest route: x_28730 + tracker x_24548 ===')
    v = list(W); v[28730] += 1000003; J.fwd2(v, 2)
    v[24548] += v[25442] - W[25442]; J.fwd2(v, 2)
    c, s, nz, av = state(v)
    print(f'  raw: out12={c} score={s} broken={nz}')
    vr, c, s, nz = repair(list(v), keep=keep_c2, tag='C2')
    print(f'  REPAIRED C2: out12={c} score={s} broken={nz}  ({time.time()-t0:.0f}s)')

    print('\n=== congruence 1, route A: x_6418 ===')
    v = list(W); v[6418] += 1000003; J.fwd2(v, 2)
    vr1, c, s, nz = repair(list(v), keep=keep_c1, tag='C1a')
    print(f'  REPAIRED C1a: out12={c} score={s} broken={nz}  ({time.time()-t0:.0f}s)')

    print('\n=== congruence 1, route B: x_7068 + tracker x_14853 ===')
    v = list(W); v[7068] += 1000003; J.fwd2(v, 2)
    v[14853] += v[1308] - W[1308]; J.fwd2(v, 2)
    vr2, c, s, nz = repair(list(v), keep=keep_c1, tag='C1b')
    print(f'  REPAIRED C1b: out12={c} score={s} broken={nz}  ({time.time()-t0:.0f}s)')
