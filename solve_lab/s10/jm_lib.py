"""Shared frame-2 helpers for the multi-parameter joint-move search (jm_*)."""
import os, sys, collections, json, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad

P = 2**256 - 2**32 - 977
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
SS = set(SEVEN)
E12 = set([2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350,
           14584, 18673, 22044, 29125])
DETACH = {7068: 22229, 28730: 22230, 29854: 35758, 31864: 35761, 642: 35762}
definer = {t: a for t, a in L.definer.items() if t not in DETACH}
atom_out = {a: o for a, o in L.atom_out.items() if a not in set(DETACH.values())}
ORDER = [t for t in ad.ORDER if t not in DETACH]
FREESET = set(u for u in range(L.NVARS) if u not in definer)
CHECKS = [a for a in range(L.NA) if a not in atom_out]
BASE = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')


def fwd2(v, rounds=2):
    for _ in range(rounds):
        for u in ORDER:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None:
                v[u] = nv
    return v


def base_state():
    v = L.load(BASE)
    w = list(v)
    fwd2(w, 8)
    return w


def resid(v):
    """the three pin residues mod p"""
    return ((v[7068] - v[2099]) % P,
            (v[14853] - v[1308]) % P,
            (v[24548] - v[25442]) % P)


def cost(v):
    """(out12, score, failing set)"""
    av = L.all_atom_values(v)
    f = set(L.failing_eqs(av))
    return len(f - E12), L.NEQ - len(f), f, av


def nz_atoms(av):
    return [a for a in range(L.NA) if av[a]]


def cone(u, dd=None):
    dd = definer if dd is None else dd
    seen = set()
    st = [u]
    while st:
        t = st.pop()
        if t in seen:
            continue
        seen.add(t)
        a = dd.get(t)
        if a is None:
            continue
        for x in L.avars[a]:
            if x != t and x not in seen:
                st.append(x)
    return seen


def jac_col(u, vm, targets):
    """forward-mode AD in frame 2: d(target var)/d(x_u) mod p for target VARS."""
    dv = collections.defaultdict(int)
    dv[u] = 1
    for t in ORDER:
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0:
            dv[t] = 0
            continue
        s = 0
        for w in L.avars[a]:
            if w == t:
                continue
            dw = dv.get(w, 0)
            if dw:
                s += ad.dpart(a, w, vm) * dw
        dv[t] = (-s % P) * pow(d, -1, P) % P
    return {t: dv.get(t, 0) % P for t in targets}


def jac_col_full(u, vm):
    dv = collections.defaultdict(int)
    dv[u] = 1
    for t in ORDER:
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0:
            dv[t] = 0
            continue
        s = 0
        for w in L.avars[a]:
            if w == t:
                continue
            dw = dv.get(w, 0)
            if dw:
                s += ad.dpart(a, w, vm) * dw
        dv[t] = (-s % P) * pow(d, -1, P) % P
    return dv
