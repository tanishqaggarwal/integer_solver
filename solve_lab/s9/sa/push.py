"""Iteratively enlarge the perturbation variable set and solve exactly over Z
for the largest subset of remaining failing equations that can be zeroed."""
import sys, itertools, json, time, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, linmodel, intsolve

_d = pickle.load(open('atoms.pkl', 'rb'))
ET = _d['eq_terms']


def apoly(a):
    return lib.roots[a] if a in lib.roots else lib.polys[a]


def conflicts(X):
    Xs = set(X)
    bad = []
    for a in range(len(lib.polys)):
        for m in apoly(a):
            if sum(1 for u in m if u in Xs) >= 2:
                bad.append((a, m))
    return bad


def cand_vars(V, F, maxextra=40):
    """Linear variables of atoms appearing in the failing equations, ranked by
    how many equations outside F they touch."""
    Fs = set(F)
    atoms = set()
    for e in F:
        for c, a in ET[e][2]:
            atoms.add(a)
    out = {}
    for a in sorted(atoms):
        for m in apoly(a):
            if len(m) != 1:
                continue
            x = m[0]
            s = set()
            for at in lib.var_atoms[x]:
                s.update(lib.atom2eq.get(at, ()))
            k = len(s - Fs)
            if k <= maxextra:
                out[x] = min(out.get(x, 10**9), k)
    return sorted(out.items(), key=lambda t: (t[1], t[0]))


def best_fix(V, X, Fs, verbose=True, cap=None):
    rows, rhs, issq, eqs = linmodel.build(X, V)
    touched = [e for e in eqs if any(rows[e])]
    must = [e for e in touched if e not in Fs]
    cand = [e for e in touched if e in Fs]
    n = len(X)
    if verbose:
        print(f'  |X|={n} touched={len(touched)} must={len(must)} cand={len(cand)}')
    top = len(cand) if cap is None else min(cap, len(cand))
    for k in range(top, 0, -1):
        for sub in itertools.combinations(cand, k):
            M = [rows[e] for e in must] + [rows[e] for e in sub]
            b = [0] * len(must) + [-rhs[e] for e in sub]
            x, info = intsolve.solve_int(M, b, n)
            if x is not None:
                return k, sub, x
        if verbose:
            print(f'    k={k}: none')
    return 0, (), None
