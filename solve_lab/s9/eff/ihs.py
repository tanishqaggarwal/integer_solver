"""Implicit-hitting-set maximiser: smallest set of equations to DROP so that the rest
of S is simultaneously zeroable over Z.  Solvability is downward closed, so the optimum
is the minimum hitting set of the minimal unsolvable cores (found lazily).
"""
import model as MD


def mhs_k(cores, k):
    """Any hitting set of size <= k, or None.  Iterative-deepening DFS."""
    if not cores:
        return set()
    if k <= 0:
        return None
    c = min(cores, key=len)
    for e in sorted(c):
        rest = [x for x in cores if e not in x]
        r = mhs_k(rest, k - 1)
        if r is not None:
            return r | {e}
    return None


def min_hitting_set(cores, cap):
    if not cores:
        return set()
    for k in range(1, cap + 1):
        r = mhs_k(cores, k)
        if r is not None:
            return r
    return None


def shrink(mod, T):
    T = set(T)
    for e in sorted(T):
        if MD.solvable(mod, T - {e}) is None:
            T -= {e}
    return frozenset(T)


def min_drop(mod, S, cap=8, maxiter=200000):
    """(dropped_set, solution) with |dropped| minimal and <= cap, else (None, None)."""
    S = set(S)
    cores = []
    for _ in range(maxiter):
        H = min_hitting_set(cores, cap)
        if H is None:
            return None, None
        T = S - H
        z = MD.solvable(mod, T)
        if z is not None:
            return H, z
        cores.append(shrink(mod, T))
    return None, None
