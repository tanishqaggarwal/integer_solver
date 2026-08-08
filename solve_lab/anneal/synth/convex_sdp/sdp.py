#!/usr/bin/env python3
"""sdp.py -- pure-numpy Shor / Lasserre level-1 SDP relaxation for QUBO/Ising.

No cvxpy, no scipy in this environment; we solve the diag=1 PSD SDP

        minimize  <C, X>     s.t.  X >= 0,  diag(X) = 1

with the MIXING METHOD (Wang, Chang & Kolter 2017): a Burer-Monteiro
coordinate descent over unit vectors  v_a in R^k.  For  k > sqrt(2N)  the
BM problem has no spurious local minima (Boumal-Voroninski-Bandeira), so the
sweep converges to the *global* SDP optimum; we run several restarts and keep
the minimum value (every feasible X is an upper bound on the SDP optimum, so
the smallest one found is the tightest estimate of it).

The SDP optimum is a LOWER bound on the true integer minimum E_min, so

        additive integrality gap  =  E_min - SDP_opt  >= 0.

gap == 0  <=>  the degree-2 SOS certificate is exact for this Hamiltonian.
"""
import numpy as np


# ---- QUBO dict -> homogenised Ising cost matrix over spins {s_0=+1, s_1..s_n} ----
def qubo_to_C(Q, n):
    """Return symmetric (n+1)x(n+1) C with  <C, s s^T> = E(x),  x_i=(1+s_i)/2,
    s_0 the +1 anchor at index 0 (so variable i lives at row/col i+1)."""
    C = np.zeros((n + 1, n + 1))
    for m, c in Q.items():
        c = float(c)
        if len(m) == 0:
            C[0, 0] += c
        elif len(m) == 1:
            i = m[0] + 1
            C[0, 0] += c / 2
            C[0, i] += c / 4
            C[i, 0] += c / 4
        else:
            i, j = m[0] + 1, m[1] + 1
            C[0, 0] += c / 4
            C[0, i] += c / 8; C[i, 0] += c / 8
            C[0, j] += c / 8; C[j, 0] += c / 8
            C[i, j] += c / 8; C[j, i] += c / 8
    return C


def sdp_min(C, k=None, restarts=24, iters=2000, seed=0, tol=1e-9):
    """Mixing-method solve of  min <C,X>, diag(X)=1, X>=0.
    Returns (value, X, V) for the best (lowest-value) restart."""
    N = C.shape[0]
    if k is None:
        k = max(3, int(np.ceil(np.sqrt(2 * N))) + 2)
    k = min(k, N)
    # adaptive budget: coordinate sweeps are Gauss-Seidel and converge fast;
    # scale work down for large N so big modmuls stay tractable.
    if N > 120:
        restarts = min(restarts, 6); iters = min(iters, 400)
    elif N > 60:
        restarts = min(restarts, 12); iters = min(iters, 700)
    rng = np.random.default_rng(seed)
    diagsum = float(np.trace(C))
    Coff = C - np.diag(np.diag(C))
    best_val, best_X, best_V = np.inf, None, None
    for r in range(restarts):
        V = rng.standard_normal((N, k))
        V /= np.linalg.norm(V, axis=1, keepdims=True)
        for _ in range(iters):
            maxchange = 0.0
            for a in range(N):
                g = Coff[a] @ V                     # sum_b Coff[a,b] v_b
                nrm = np.linalg.norm(g)
                if nrm < 1e-13:
                    nv = rng.standard_normal(k); nv /= np.linalg.norm(nv)
                else:
                    nv = -g / nrm                   # minimise <v_a, g>
                d = nv - V[a]
                cc = d @ d
                if cc > maxchange:
                    maxchange = cc
                V[a] = nv
            if maxchange < tol:
                break
        X = V @ V.T
        val = float(np.sum(Coff * X)) + diagsum
        if val < best_val:
            best_val, best_X, best_V = val, X, V.copy()
    return best_val, best_X, best_V


def rank_of(X, rel=1e-6):
    """numerical rank via eigenvalues of the symmetric PSD X."""
    w = np.linalg.eigvalsh((X + X.T) / 2)
    w = np.clip(w, 0, None)
    thr = rel * w.max() if w.max() > 0 else 0
    return int((w > thr).sum()), w[::-1]


def brute_min(Q, n):
    """exact ground energy and full ground-state list by 2^n enumeration."""
    def E(x):
        e = 0
        for m, c in Q.items():
            if not m: e += c
            elif len(m) == 1: e += c * x[m[0]]
            else: e += c * x[m[0]] * x[m[1]]
        return e
    best = None; states = []
    for code in range(1 << n):
        x = [(code >> i) & 1 for i in range(n)]
        e = E(x)
        if best is None or e < best:
            best = e; states = [x]
        elif e == best:
            states.append(x)
    return best, states
