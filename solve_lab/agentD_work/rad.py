"""Reverse-mode AD over Z/qZ (q a large random prime, NOT p) for the gate DAG.

Gives, for a target atom c, the derivative d c / d u for every variable u, in one
backward pass.  Used only to shortlist candidate knobs; every proposed move is
then verified by exact integer probing.
"""
import collections, random, sys
import dlib as L

Q = (1 << 61) - 1   # Mersenne prime 2305843009213693951


def datom(a, v, q=Q):
    """Partial derivatives of atom a wrt each of its variables, mod q."""
    d = collections.defaultdict(int)
    for m, c in L.polys[a].items():
        n = len(m)
        if n == 0:
            continue
        for j in range(n):
            t = c % q
            for k in range(n):
                if k != j:
                    t = t * (v[m[k]] % q) % q
            d[m[j]] = (d[m[j]] + t) % q
    return d


_rtopo = None


def rtopo():
    global _rtopo
    if _rtopo is None:
        _rtopo = list(reversed(L.topo))
    return _rtopo


def grad(c, v, q=Q):
    """Return dict var -> d(atom c)/d(var) mod q, propagated to free inputs."""
    adj = collections.defaultdict(int)
    for u, dv in datom(c, v, q).items():
        adj[u] = (adj[u] + dv) % q
    for t in rtopo():
        at = adj.get(t)
        if not at:
            continue
        a = L.definer[t]
        d = datom(a, v, q)
        ct = d.get(t, 0)
        if ct == 0:
            continue
        inv = pow(ct, q - 2, q)
        for w, dw in d.items():
            if w == t or dw == 0:
                continue
            adj[w] = (adj[w] - at * dw % q * inv) % q
    return adj


def free_knobs(c, v, q=Q):
    g = grad(c, v, q)
    return {u: d for u, d in g.items() if d and u in L.freeset}
