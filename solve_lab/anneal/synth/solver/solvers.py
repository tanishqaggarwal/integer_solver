#!/usr/bin/env python3
"""solvers.py -- a battery of classical Ising solvers, all on the numpy Ising model
from model.py.  Every solver returns (best_energy, best_spin_state) where best_spin
is a {-1,+1} numpy array; the QUBO energy is recovered as ising.energy(best_spin),
which is 0 exactly on a ground state.

Solvers:
  sa        simulated annealing (geometric beta ramp, incremental grad)     -- baseline
  pt        parallel tempering / replica exchange across a beta ladder
  tabu      steepest 1-opt descent + tabu list + 2-opt escape (qbsolv inner loop)
  sb        simulated bifurcation (Toshiba), ballistic (bSB) or discrete (dSB)

`clamp` is an iterable of variable indices held fixed at their value in `x0`
(spins frozen).  All solvers respect it.
"""
import numpy as np


def _clamp_arrays(ising, clamp, x0):
    n = ising.n
    fixed = np.zeros(n, dtype=bool)
    s = np.empty(n)
    if x0 is not None:
        s[:] = 2 * np.asarray(x0, dtype=np.float64) - 1
    else:
        s[:] = 1.0
    if clamp:
        idx = np.fromiter(clamp, dtype=np.int64)
        fixed[idx] = True
    free = np.where(~fixed)[0]
    return s, fixed, free


# ---------------------------------------------------------------- SA
def sa(ising, sweeps=2000, beta0=0.05, beta1=30.0, seed=0, x0=None, clamp=(), rng=None):
    rng = rng or np.random.default_rng(seed)
    n = ising.n
    s, fixed, free = _clamp_arrays(ising, clamp, x0)
    if x0 is None:
        s[free] = rng.integers(0, 2, size=len(free)) * 2 - 1
    indptr, indices, data, h = ising.indptr, ising.indices, ising.data, ising.h
    grad = ising.grad(s)
    e = ising.energy(s)
    best, bs = e, s.copy()
    betas = beta0 * (beta1 / beta0) ** (np.arange(sweeps) / max(1, sweeps - 1))
    for sw in range(sweeps):
        beta = betas[sw]
        order = free.copy(); rng.shuffle(order)
        u = rng.random(len(order))
        for idx, ui in zip(order, u):
            gi = grad[idx]
            si = s[idx]
            de = -2.0 * si * gi
            if de <= 0 or ui < np.exp(-beta * de):
                ds = -2.0 * si
                s[idx] = -si
                e += de
                a, b = indptr[idx], indptr[idx + 1]
                grad[indices[a:b]] += ds * data[a:b]
                if e < best:
                    best = e; bs = s.copy()
        if best <= 1e-9:
            break
    return best, bs


# ---------------------------------------------------------------- Parallel tempering
def pt(ising, sweeps=2000, R=12, beta0=0.05, beta1=30.0, seed=0, x0=None, clamp=(),
       rng=None, swap_every=1):
    rng = rng or np.random.default_rng(seed)
    n = ising.n
    betas = beta0 * (beta1 / beta0) ** (np.arange(R) / max(1, R - 1))
    s0, fixed, free = _clamp_arrays(ising, clamp, x0)
    S = np.tile(s0, (R, 1))
    if x0 is None:
        S[:, free] = rng.integers(0, 2, size=(R, len(free))) * 2 - 1
    else:
        # diversify replicas on the free set
        for r in range(R):
            S[r, free] = rng.integers(0, 2, size=len(free)) * 2 - 1
        S[:, list(clamp)] = s0[list(clamp)] if clamp else S[:, :0]
    indptr, indices, data, h = ising.indptr, ising.indices, ising.data, ising.h
    grads = [ising.grad(S[r]) for r in range(R)]
    energies = np.array([ising.energy(S[r]) for r in range(R)])
    best = energies.min(); bs = S[int(energies.argmin())].copy()
    for sw in range(sweeps):
        for r in range(R):
            beta = betas[r]
            s = S[r]; grad = grads[r]
            order = free.copy(); rng.shuffle(order)
            u = rng.random(len(order))
            e = energies[r]
            for idx, ui in zip(order, u):
                si = s[idx]; gi = grad[idx]
                de = -2.0 * si * gi
                if de <= 0 or ui < np.exp(-beta * de):
                    ds = -2.0 * si
                    s[idx] = -si; e += de
                    a, b = indptr[idx], indptr[idx + 1]
                    grad[indices[a:b]] += ds * data[a:b]
            energies[r] = e
            if e < best:
                best = e; bs = s.copy()
        if sw % swap_every == 0:
            # replica exchange on adjacent temperatures
            start = sw % 2
            for r in range(start, R - 1, 2):
                d = (betas[r] - betas[r + 1]) * (energies[r] - energies[r + 1])
                if d >= 0 or rng.random() < np.exp(d):
                    S[[r, r + 1]] = S[[r + 1, r]]
                    grads[r], grads[r + 1] = grads[r + 1], grads[r]
                    energies[[r, r + 1]] = energies[[r + 1, r]]
        if best <= 1e-9:
            break
    return best, bs


# ---------------------------------------------------------------- Tabu (1-opt + 2-opt)
def tabu(ising, iters=20000, tenure=None, seed=0, x0=None, clamp=(), rng=None,
         restarts=1):
    rng = rng or np.random.default_rng(seed)
    n = ising.n
    best_global = None; bs_global = None
    tenure = tenure or max(10, n // 100)
    indptr, indices, data, h = ising.indptr, ising.indices, ising.data, ising.h
    for rs in range(restarts):
        s, fixed, free = _clamp_arrays(ising, clamp, x0)
        if x0 is None or rs > 0:
            s[free] = rng.integers(0, 2, size=len(free)) * 2 - 1
        grad = ising.grad(s)
        e = ising.energy(s)
        best, bs = e, s.copy()
        tabu_until = np.zeros(n, dtype=np.int64)
        freemask = ~fixed
        it = 0
        stall = 0
        while it < iters:
            # delta for every free var: dE_i = -2 s_i grad_i
            delta = -2.0 * s * grad
            # forbid clamped
            delta_masked = np.where(freemask, delta, np.inf)
            # tabu: raise cost unless it would beat global best (aspiration)
            tmask = (tabu_until > it) & (e + delta_masked >= best - 1e-9)
            delta_masked = np.where(tmask, np.inf, delta_masked)
            idx = int(np.argmin(delta_masked))
            if not np.isfinite(delta_masked[idx]):
                # everything tabu; do a random kick (2-opt style double flip)
                cand = rng.choice(free, size=2, replace=False)
                for idx in cand:
                    si = s[idx]; ds = -2.0 * si
                    s[idx] = -si; e += -2.0 * si * grad[idx]
                    a, b = indptr[idx], indptr[idx + 1]
                    grad[indices[a:b]] += ds * data[a:b]
                it += 1
                continue
            si = s[idx]; ds = -2.0 * si
            s[idx] = -si; e += delta[idx]
            a, b = indptr[idx], indptr[idx + 1]
            grad[indices[a:b]] += ds * data[a:b]
            tabu_until[idx] = it + tenure
            if e < best - 1e-9:
                best = e; bs = s.copy(); stall = 0
            else:
                stall += 1
            it += 1
            if best <= 1e-9:
                break
            if stall > 4 * tenure:
                # diversification kick: flip a handful of random free spins
                cand = rng.choice(free, size=min(len(free), 5), replace=False)
                for k in cand:
                    si = s[k]; ds = -2.0 * si
                    s[k] = -si; e += -2.0 * si * grad[k]
                    a, b = indptr[k], indptr[k + 1]
                    grad[indices[a:b]] += ds * data[a:b]
                stall = 0
        if best_global is None or best < best_global:
            best_global = best; bs_global = bs.copy()
        if best_global <= 1e-9:
            break
    return best_global, bs_global


# ---------------------------------------------------------------- Simulated bifurcation
def sb(ising, steps=1000, dt=0.5, seed=0, x0=None, clamp=(), rng=None,
       kind='dSB', a0=1.0, c0=None, restarts=1):
    """Toshiba simulated bifurcation.  kind='dSB' (discrete, sign(x) in coupling)
    or 'bSB' (ballistic, x in coupling).  Fully vectorised."""
    rng = rng or np.random.default_rng(seed)
    n = ising.n
    if c0 is None:
        c0 = 0.5 / (ising.Jstd * np.sqrt(max(1, n)))
    s0, fixed, free = _clamp_arrays(ising, clamp, x0)
    clamp_spin = s0.copy()
    best_global = None; bs_global = None
    for rs in range(restarts):
        x = 0.1 * (rng.random(n) - 0.5)
        y = 0.1 * (rng.random(n) - 0.5)
        # clamp fixed vars to their spin
        x[fixed] = clamp_spin[fixed]
        y[fixed] = 0.0
        for t in range(steps):
            a_t = a0 * (t + 1) / steps
            if kind == 'dSB':
                phi = np.sign(x)
                phi[phi == 0] = 1.0
                f = ising.h + ising.Amatvec(phi)
            else:  # bSB
                f = ising.h + ising.Amatvec(x)
            y += dt * (-(a0 - a_t) * x - c0 * f)
            x += dt * a0 * y
            # inelastic walls
            over = np.abs(x) > 1.0
            x[over] = np.sign(x[over])
            y[over] = 0.0
            # re-pin clamped
            x[fixed] = clamp_spin[fixed]
            y[fixed] = 0.0
        spin = np.sign(x); spin[spin == 0] = 1.0
        spin[fixed] = clamp_spin[fixed]
        e = ising.energy(spin)
        # local 1-opt polish (SB output is often 1-2 flips from a min)
        e, spin = _greedy_polish(ising, spin, fixed)
        if best_global is None or e < best_global:
            best_global = e; bs_global = spin.copy()
        if best_global <= 1e-9:
            break
    return best_global, bs_global


def _greedy_polish(ising, s, fixed, max_pass=50):
    s = s.copy()
    grad = ising.grad(s)
    e = ising.energy(s)
    indptr, indices, data = ising.indptr, ising.indices, ising.data
    freemask = ~fixed
    for _ in range(max_pass):
        delta = -2.0 * s * grad
        delta = np.where(freemask, delta, np.inf)
        idx = int(np.argmin(delta))
        if delta[idx] >= -1e-12:
            break
        si = s[idx]; ds = -2.0 * si
        s[idx] = -si; e += delta[idx]
        a, b = indptr[idx], indptr[idx + 1]
        grad[indices[a:b]] += ds * data[a:b]
    return e, s
