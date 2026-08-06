"""Handle-set construction helpers (importable; no side effects)."""
import sys, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib

ET = pickle.load(open('atoms.pkl', 'rb'))['eq_terms']


def apoly(a):
    return lib.roots[a] if a in lib.roots else lib.polys[a]


def ok_to_add(X, x):
    """No monomial of any atom may contain x together with another member of X
    (or x twice) -- otherwise the linear model is not exact."""
    Xs = set(X)
    for a in lib.var_atoms.get(x, ()):
        for m in apoly(a):
            if x not in m:
                continue
            if m.count(x) >= 2:
                return False
            if any(u in Xs for u in m if u != x):
                return False
    return True


def eff_foot(x, V):
    """Equations actually moved by perturbing x at the state V."""
    E = set()
    for a in lib.var_atoms.get(x, ()):
        Pp = apoly(a)
        d = 0
        for m, c in Pp.items():
            k = m.count(x)
            if k == 0:
                continue
            t = c * k
            rem = list(m)
            rem.remove(x)
            for u in rem:
                t *= V[u]
            d += t
        if d:
            E.update(lib.atom2eq.get(a, ()))
    return E


def vars_of_eqs(eqs):
    out = set()
    for e in eqs:
        for c, a in ET[e][2]:
            for m in apoly(a):
                out.update(m)
    return out


def build_handles(V, Fs, nvars, seedv=()):
    X = []
    for x in seedv:
        if ok_to_add(X, x):
            X.append(x)
    L1 = vars_of_eqs(sorted(Fs))
    seedE = set(Fs)
    for x in L1:
        seedE |= eff_foot(x, V)
    cands = sorted(L1 | vars_of_eqs(sorted(seedE)))
    scored = []
    for x in cands:
        if x in X:
            continue
        E = eff_foot(x, V)
        if not E:
            continue
        scored.append((len(E - Fs), -len(E & Fs), x))
    scored.sort()
    for _, _, x in scored:
        if len(X) >= nvars:
            break
        if ok_to_add(X, x):
            X.append(x)
    return X
