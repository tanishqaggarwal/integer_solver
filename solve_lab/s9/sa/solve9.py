"""Exact rational solve: perturb a chosen variable set to kill as many of the
remaining failing equations as possible while keeping every other equation zero."""
import sys, pickle, itertools, json
from fractions import Fraction
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib, linmodel

V = list(lib.V0)
V[28730] = V[4432] - V[19964]          # the 10-failing state
F = lib.true_fails(V)
print('start state failing:', len(F), F)
Fs = set(F)


def rref(M, ncol):
    """M: list of rows (list of Fraction, length ncol+1 incl rhs). In-place RREF.
    Returns list of pivot columns."""
    piv = []
    r = 0
    for c in range(ncol):
        pr = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                pr = i; break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    return piv


def solve_system(rows, rhs, X, must_zero, want_zero):
    """Find rational D with row.D = 0 for e in must_zero, row.D = -rhs for e in want_zero.
    Returns D (list of Fraction) or None."""
    n = len(X)
    M = []
    for e in must_zero:
        M.append([Fraction(c) for c in rows[e]] + [Fraction(0)])
    for e in want_zero:
        M.append([Fraction(c) for c in rows[e]] + [Fraction(-rhs[e])])
    piv = rref(M, n)
    # inconsistent?
    for row in M:
        if all(row[c] == 0 for c in range(n)) and row[n] != 0:
            return None
    D = [Fraction(0)] * n
    for i, c in enumerate(piv):
        D[c] = M[i][n]
    return D


def attempt(X, maxfix=None, verbose=True):
    rows, rhs, issq, eqs = linmodel.build(X, V)
    touched = [e for e in eqs if any(rows[e])]
    must = [e for e in touched if e not in Fs]
    cand = [e for e in touched if e in Fs]
    if verbose:
        print(f'  vars={len(X)} touched eqs={len(touched)} must-keep={len(must)} fixable-candidates={len(cand)}')
    best = None
    top = len(cand) if maxfix is None else min(maxfix, len(cand))
    for k in range(top, 0, -1):
        for sub in itertools.combinations(cand, k):
            D = solve_system(rows, rhs, X, must, list(sub))
            if D is not None:
                return sub, D, rows, rhs, must, cand
        if k <= top - 3 and best is None:
            pass
    return None, None, rows, rhs, must, cand
