"""S10 shared tools: effective-linear atom solving + equation-space scored repair."""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
NVARS = L.NVARS


def lin_parts(a, u, v):
    """Return (c, rest) with atom_a(v) == c*v[u] + rest, valid when x_u appears at
    most to the first power in every monomial of atom a.  Returns None if x_u
    appears squared."""
    c = 0
    rest = 0
    for m, cc in L.polys[a].items():
        k = m.count(u)
        if k > 1:
            return None
        t = cc
        if k == 0:
            for w in m:
                t *= v[w]
            rest += t
        else:
            for w in m:
                if w != u:
                    t *= v[w]
            c += t
    return c, rest


def solve_lin(a, u, v):
    """Integer value for x_u making atom a vanish, or None."""
    r = lin_parts(a, u, v)
    if r is None:
        return None
    c, rest = r
    if c == 0 or rest % c:
        return None
    return -rest // c


def atom_vars(a):
    return L.avars[a]


def measure(v, tag=None, show=0):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    if tag:
        print(f'[{tag}] nz_atoms={len(nz)} failing={len(fail)} score={L.NEQ - len(fail)}')
        for a in nz[:show]:
            print(f'    a{a:<6} neq={len(L.atom2eq.get(a, {})):<3} {L.atom_src[a][:100]}')
    return av, nz, fail


def try_fix(v, a, block=(), verbose=False):
    """Try every variable of atom `a`; return list of (score, var, newval) sorted best-first."""
    out = []
    for u in sorted(atom_vars(a)):
        nv = solve_lin(a, u, v)
        if nv is None or nv == v[u]:
            continue
        w = list(v)
        try:
            L.ripple(w, {u: nv}, block=block)
        except Exception:
            continue
        av = L.all_atom_values(w)
        sc = L.NEQ - len(L.failing_eqs(av))
        out.append((sc, u, nv))
        if verbose:
            print(f'      fix a{a} via x_{u}: score {sc}')
    out.sort(key=lambda t: -t[0])
    return out


def save(v, path):
    json.dump({f'x_{i}': v[i] for i in range(NVARS) if v[i] != 0}, open(path, 'w'))
