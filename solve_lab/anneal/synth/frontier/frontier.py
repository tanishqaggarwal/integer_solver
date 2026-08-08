#!/usr/bin/env python3
"""frontier.py -- measure F, the arithmetic-annealing frontier: the largest field
size s (bits) whose modular multiply a*b==c (mod p) a real annealer drives to the
ground state (QUBO energy exactly 0) with meaningful probability.

Success probability = fraction of INDEPENDENT restarts reaching E=0.
Every cell is time-capped and checkpointed to JSON as it completes.
"""
import sys, os, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'synth', 'solver'))
sys.path.insert(0, ROOT)
import numpy as np
import solvers as S
import fbuild as FB


# ----------------------------------------------------------- effort tiers
EFFORT = {
    'low':     dict(sa=2000,   pt=(1000, 8),   tabu=8000,   sb=1500),
    'mid':     dict(sa=20000,  pt=(6000, 12),  tabu=40000,  sb=8000),
    'high':    dict(sa=80000,  pt=(20000, 16), tabu=150000, sb=30000),
    'extreme': dict(sa=250000, pt=(60000, 24), tabu=600000, sb=120000),
}


def _spin_to_bits(spin):
    return ((np.asarray(spin) + 1) // 2).astype(int).tolist()


def run_one(ising, solver, budget, seed, x0=None, clamp=()):
    """one independent restart -> best energy reached."""
    if solver == 'sa':
        e, _ = S.sa(ising, sweeps=budget['sa'], seed=seed, x0=x0, clamp=clamp)
    elif solver == 'pt':
        sw, R = budget['pt']
        e, _ = S.pt(ising, sweeps=sw, R=R, seed=seed, x0=x0, clamp=clamp)
    elif solver == 'tabu':
        e, _ = S.tabu(ising, iters=budget['tabu'], seed=seed, x0=x0, clamp=clamp)
    elif solver == 'sb':
        e, _ = S.sb(ising, steps=budget['sb'], seed=seed, x0=x0, clamp=clamp)
    elif solver == 'sb_tabu':
        # simulated bifurcation, then tabu polish from its output state
        e1, sp = S.sb(ising, steps=budget['sb'], seed=seed, x0=x0, clamp=clamp)
        if e1 == 0:
            return 0.0
        e, _ = S.tabu(ising, iters=budget['tabu'], seed=seed,
                      x0=_spin_to_bits(sp), clamp=clamp)
        e = min(e, e1)
    elif solver == 'pt_wide':
        # parallel tempering with a doubled replica ladder (barrier crosser)
        sw, R = budget['pt']
        e, _ = S.pt(ising, sweeps=sw, R=2 * R, seed=seed, x0=x0, clamp=clamp)
    else:
        raise ValueError(solver)
    return e


def success_rate(ising, solver, budget, n_max, time_cap, seed0=0,
                 x0=None, clamp=()):
    """run up to n_max independent restarts or until time_cap seconds; return
    (hits, n_done, best_energy, elapsed)."""
    hits = 0
    best = 1e18
    t0 = time.time()
    n = 0
    for r in range(n_max):
        e = run_one(ising, solver, budget, seed0 + r, x0=x0, clamp=clamp)
        best = min(best, e)
        hits += (e == 0)
        n += 1
        if time.time() - t0 > time_cap:
            break
    return hits, n, float(best), time.time() - t0


# ----------------------------------------------------------- driver
def sweep(sizes, solvers, efforts, encoder='baseline', mode='wallace',
          W_and=None, n_caps=None, time_cap=45.0, clamp_which=None,
          seed0=1000, out_json=None, squeeze_kw=None, label=''):
    """Full solver x s x effort table.  n_caps: {solver: n_max}."""
    n_caps = n_caps or dict(sa=40, pt=24, pt_wide=16, tabu=400, sb=600,
                            sb_tabu=300)
    results = {}
    for s in sizes:
        if encoder == 'baseline':
            mm = FB.baseline_modmul(s, mode=mode, W_and=W_and)
        else:
            mm = FB.squeeze_modmul(s, **(squeeze_kw or {}))
        x0, clamp = (None, ())
        if clamp_which:
            x0, clamp = FB.clamp_operand(mm, which=clamp_which,
                                         rng=np.random.default_rng(7))
        for eff in efforts:
            budget = EFFORT[eff]
            for solver in solvers:
                key = f"{s}|{eff}|{solver}"
                hits, n, best, el = success_rate(
                    mm['ising'], solver, budget, n_caps.get(solver, 100),
                    time_cap, seed0=seed0, x0=x0, clamp=clamp)
                hits = int(hits); n = int(n)
                rate = hits / n if n else 0.0
                results[key] = dict(s=int(s), effort=eff, solver=solver,
                                    n_vars=int(mm['Q'].n),
                                    hits=hits, n=n, rate=float(rate),
                                    best=float(best), elapsed=round(el, 1))
                print(f"[{label}] s={s:2d} n={mm['Q'].n:5d} {eff:7s} {solver:8s} "
                      f"{hits:4d}/{n:<4d} rate={rate:6.3f} best={best:5.0f} ({el:4.1f}s)")
                if out_json:
                    json.dump(results, open(out_json, 'w'), indent=1)
    return results


if __name__ == '__main__':
    pass
