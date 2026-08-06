"""Exact linear model of equation residuals under direct perturbation of a chosen
variable set (gate outputs included -- i.e. we are allowed to BREAK gates).

Each atom is degree <= 2.  If no chosen variable appears quadratically and no two
chosen variables multiply each other inside one atom, then
    atom_a(v + D) = atom_a(v) + sum_x d[a][x] * D[x]      (EXACT, not first order)
with d[a][x] = dP_a/dx evaluated at v.
Equation residual (non-square) = m_e * sum_a c[e][a] * atom_a.
Square equation residual = (root_e)^2, modelled through root_e.
"""
import sys, pickle
from fractions import Fraction
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib

_d = pickle.load(open('atoms.pkl', 'rb'))
EQ_TERMS = _d['eq_terms']
ATOM_SRC = _d['atom_src']


def datom(a, x, v, use_root=False):
    """dP/dx at v for atom a (or its root)."""
    P = lib.roots[a] if (use_root and a in lib.roots) else lib.polys[a]
    s = 0
    for m, c in P.items():
        k = m.count(x)
        if k == 0:
            continue
        t = c * k
        rem = list(m)
        rem.remove(x)
        for u in rem:
            t *= v[u]
        s += t
    return s


def nonlinear_conflict(X, v):
    """Return list of (atom, monomial) where the linear model would be wrong."""
    Xs = set(X)
    bad = []
    for a in range(len(lib.polys)):
        for m in lib.polys[a]:
            k = sum(1 for u in m if u in Xs)
            if k >= 2:
                bad.append((a, m))
    return bad


def build(X, v):
    """Return (rows, rhs, eqlist) where row[e][i] is the coefficient of D[X[i]]
    in the residual of equation e, and rhs[e] is the current residual."""
    Xs = {x: i for i, x in enumerate(X)}
    # atoms touched
    atoms = set()
    for x in X:
        atoms.update(lib.var_atoms[x])
    # equations touched
    eqs = set()
    for a in atoms:
        eqs.update(lib.atom2eq.get(a, ()))
    eqs = sorted(eqs)
    # atom derivative table
    D = {}
    for a in atoms:
        D[a] = {x: datom(a, x, v) for x in X if datom(a, x, v) != 0}
    rows = {}
    rhs = {}
    issq = {}
    for e in eqs:
        m, sq, tl = EQ_TERMS[e]
        r = [0] * len(X)
        cur = 0
        for c, a in tl:
            if sq:
                # single square atom: model the root
                val = lib.evalpoly(lib.roots[a] if a in lib.roots else lib.polys[a], v)
                cur += c * val
                for x in X:
                    dd = datom(a, x, v, use_root=True)
                    if dd:
                        r[Xs[x]] += c * dd
            else:
                cur += c * lib.evalpoly(lib.polys[a], v)
                for x, dd in D.get(a, {}).items():
                    r[Xs[x]] += c * dd
        rows[e] = r
        rhs[e] = cur          # residual (pre outer multiplier / pre square)
        issq[e] = sq
    return rows, rhs, issq, eqs
