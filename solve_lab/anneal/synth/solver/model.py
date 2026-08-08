#!/usr/bin/env python3
"""model.py -- build the comb / modmul QUBOs from synthetic planted-key instances
and expose them in a numpy Ising form the solvers can chew on fast.

Nothing here modifies the existing encoder; it only calls ladder.build_win /
Ladder.mul_eq and reshapes the resulting QUBO dict.

QUBO->Ising:  x_i in {0,1},  s_i = 2 x_i - 1 in {-1,+1}.
  E(x) = C0 + sum_i a_i x_i + sum_{i<j} b_ij x_i x_j
       = const + h . s + 0.5 s^T A s          (A symmetric, diag 0)
with  A_ij = b_ij/4,  h_i = a_i/2 + sum_j b_ij/4,
      const = C0 + sum_i a_i/2 + sum_{i<j} b_ij/4.
"""
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import numpy as np
from ladder import build_win, Ladder
from qubo import QB


# --------------------------------------------------------------- Ising wrapper
class Ising:
    """Symmetric-CSR Ising model with optional clamped spins."""
    def __init__(self, n, h, rows, cols, w, const):
        self.n = n
        self.h = np.asarray(h, dtype=np.float64)
        # build symmetric CSR (both (i,j) and (j,i))
        I = np.concatenate([rows, cols])
        J = np.concatenate([cols, rows])
        W = np.concatenate([w, w]).astype(np.float64)
        order = np.argsort(I, kind='stable')
        I, J, W = I[order], J[order], W[order]
        self.indptr = np.zeros(n + 1, dtype=np.int64)
        np.add.at(self.indptr, I + 1, 1)
        np.cumsum(self.indptr, out=self.indptr)
        self.indices = J.astype(np.int64)
        self.data = W
        self.const = float(const)
        # squared coupling scale for SB normalisation
        self.Jstd = float(np.sqrt((W ** 2).sum() / max(1, n))) if len(W) else 1.0

    def Amatvec(self, s):
        """returns A @ s   (sum_j A_ij s_j)."""
        out = np.zeros(self.n)
        # segment sum via reduceat is awkward with empty rows; use bincount
        contrib = self.data * s[self.indices]
        np.add.at(out, np.repeat(np.arange(self.n), np.diff(self.indptr)), contrib)
        return out

    def grad(self, s):
        return self.h + self.Amatvec(s)

    def energy(self, s):
        return self.const + float(self.h @ s) + 0.5 * float(s @ self.Amatvec(s))


def qubo_to_ising(Q, n):
    C0 = 0.0
    h = np.zeros(n)
    rows, cols, w = [], [], []
    for m, c in Q.items():
        if not m:
            C0 += c
        elif len(m) == 1:
            h[m[0]] += c
        else:
            i, j = m
            rows.append(i); cols.append(j); w.append(c)
    rows = np.array(rows, dtype=np.int64) if rows else np.zeros(0, np.int64)
    cols = np.array(cols, dtype=np.int64) if cols else np.zeros(0, np.int64)
    w = np.array(w, dtype=np.float64) if len(w) else np.zeros(0)
    # accumulate into h and const
    const = C0 + 0.5 * h.sum() + 0.25 * w.sum()
    hi = 0.5 * h.copy()
    # h_i += sum_j b_ij/4  (over both endpoints)
    np.add.at(hi, rows, 0.25 * w)
    np.add.at(hi, cols, 0.25 * w)
    A_w = 0.25 * w
    return Ising(n, hi, rows, cols, A_w, const)


# --------------------------------------------------------------- comb builder
def build_comb(inst, mu, w=2, mode='wallace', neq=False, W_and=None, secret=None):
    """Build the windowed-comb QUBO for a mu-bit residual dlog on inst's curve.

    Plants `secret` (default inst.k mod 2^mu) as the answer; returns everything a
    solver needs: the QB, the Ising, the one-hot answer vars, the planted digits,
    and the exact ground state (E=0 witness)."""
    c, G, n = inst.curve, inst.G, inst.n
    if secret is None:
        secret = inst.k & ((1 << mu) - 1)
    M = (mu + w - 1) // w
    table = [[c.mul(((t + 1) << (w * j)) % n, G) for t in range(1 << w)] for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))
    Tp = c.add(c.mul(secret % n, G), c.mul(off % n, G))
    if W_and is not None:
        # monkey-free: build then re-finalize with a chosen AND weight
        L = Ladder(c.p, mode=mode)
        L.qb.W_and = W_and
        # rebuild via build_win with W_and set: build_win makes its own Ladder,
        # so instead we post-process (see build_win path below)
    L, U = build_win(c.p, c.B, table, Tp, w, mode=mode, neq=neq)
    Q = L.qb
    if W_and is not None:
        # recompute Q dict with a different AND weight, without touching the file
        Q.W_and = W_and
        Q.finalize()
    digits = [(secret >> (w * j)) % (1 << w) for j in range(M)]
    # ground-state witness
    wv0 = {f"_u{j}": digits[j] for j in range(M)}
    inp = {}
    for j in range(M):
        for t in range(1 << w):
            inp[U[j][t]] = 1 if t == digits[j] else 0
    xstar, _ = Q.witness(inp, wv0)
    assert Q.energy(xstar) == 0, "planted secret is not a zero-energy state"
    answer_vars = [U[j][t] for j in range(M) for t in range(1 << w)]
    ising = qubo_to_ising(Q.Q, Q.n)
    # sanity: ising energy of ground state == 0
    sstar = 2 * np.array(xstar, dtype=np.float64) - 1
    e_is = ising.energy(sstar)
    assert abs(e_is) < 1e-6, f"ising ground energy {e_is} != 0"
    return dict(Q=Q, U=U, M=M, w=w, secret=secret, digits=digits, xstar=xstar,
                answer_vars=answer_vars, ising=ising, inst=inst, mu=mu)


def onehot_start(model, digits=None, rng=None):
    """a full x0 with the one-hot answer set to `digits` (default planted) and
    every other bit random. Returns (x0 list, clamp set of answer vars)."""
    Q, U, M, w = model['Q'], model['U'], model['M'], model['w']
    digits = model['digits'] if digits is None else digits
    x0 = [0] * Q.n
    if rng is not None:
        rbit = (lambda: int(rng.integers(0, 2))) if hasattr(rng, 'integers') else (lambda: rng.randrange(2))
        for i in range(Q.n):
            x0[i] = rbit()
    for j in range(M):
        for t in range(1 << w):
            x0[U[j][t]] = 1 if t == digits[j] else 0
    clamp = set(model['answer_vars'])
    return x0, clamp


# --------------------------------------------------------------- modmul unit
def build_modmul(s, mode='wallace', seed=3, W_and=None):
    """One modular multiplication a*b == c (mod p) as a QUBO (the atom, cf unit_probe).
    Returns QB, Ising, the input-word var lists A,B, and the planted values."""
    import random
    rnd = random.Random(seed)
    p = (1 << (s - 1)) + 2 * rnd.randrange(1 << (s - 3)) + 1
    while True:
        if all(p % q for q in range(3, int(p ** .5) + 1, 2)):
            break
        p += 2
    a, b = rnd.randrange(2, p), rnd.randrange(2, p)
    cc = a * b % p
    L = Ladder(p, mode=mode)
    Q = L.qb
    if W_and is not None:
        Q.W_and = W_and
    A = Q.word("a", s, lambda wv: a)
    Bw = Q.word("b", s, lambda wv: b)
    L.mul_eq("mm", A, Bw, "a", "b", [], cc)
    Q.finalize()
    x, _ = Q.witness({}, {})
    assert Q.energy(x) == 0
    ising = qubo_to_ising(Q.Q, Q.n)
    return dict(Q=Q, ising=ising, A=A, B=Bw, p=p, a=a, b=b, c=cc, xstar=x, s=s)
