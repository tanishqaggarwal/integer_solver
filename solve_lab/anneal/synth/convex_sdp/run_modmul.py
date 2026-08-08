#!/usr/bin/env python3
"""run_modmul.py -- SDP relaxation of one modular multiply a*b == c (mod p).

For s = 4,6,8:
  * enumerate ALL ground states by sweeping the 2^(2s) free (a,b) input words
    and replaying the encoder's witness to fill every ancilla (soundness gate);
  * solve the Shor SDP and measure the additive integrality gap vs E_min=0;
  * detect SDP-persistent variables (|X[0,i]|~1) and VERIFY each against the
    enumerated ground states before counting it removable.
"""
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import numpy as np
from synth.solver.model import build_modmul
from sdp import qubo_to_C, sdp_min, rank_of


def enumerate_ground_states(model):
    """Return (list of full x assignments with E==0, per-var pinned info).
    Overrides the planted 'a'/'b' word fns so witness replays an arbitrary
    input pair; keeps only pairs that replay cleanly to zero energy."""
    Q, s = model['Q'], model['s']
    box = {'a': model['a'], 'b': model['b']}
    # monkeypatch the two input-word fns on THIS instance's trace (not the file)
    for i, op in enumerate(Q.trace):
        if op[0] == 'word' and op[1] == 'a':
            Q.trace[i] = ('word', 'a', op[2], lambda wv: box['a'])
        elif op[0] == 'word' and op[1] == 'b':
            Q.trace[i] = ('word', 'b', op[2], lambda wv: box['b'])
    states = []
    for a in range(1 << s):
        for b in range(1 << s):
            box['a'], box['b'] = a, b
            try:
                x, _ = Q.witness({}, {})
            except Exception:
                continue
            if Q.energy(x) == 0:
                states.append(x)
    arr = np.array(states)
    pinned = {i: int(arr[0, i]) for i in range(Q.n) if len(set(arr[:, i])) == 1}
    return states, pinned


def run(s):
    model = build_modmul(s, mode='wallace', seed=3)
    Q = model['Q']; n = Q.n
    states, pinned = enumerate_ground_states(model)
    C = qubo_to_C(Q.Q, n)
    val, X, V = sdp_min(C, restarts=20, iters=1500, seed=1)
    gap = 0.0 - val                       # E_min == 0 for a feasible instance
    rk, eig = rank_of(X)
    # SDP-persistent candidates: |correlation with anchor| ~ 1
    corr = X[0, 1:]                        # X[0, i+1] = <s_0, s_i>
    sdp_cand = [i for i in range(n) if abs(corr[i]) > 0.999]
    # verify each candidate against enumerated ground states
    verified = [i for i in sdp_cand if i in pinned]
    print(f"=== modmul s={s}  p={model['p']}  n_vars={n}  couplers={sum(1 for k in Q.Q if len(k)==2)}")
    print(f"    ground states (E=0): {len(states)}")
    print(f"    SDP_opt={val:+.5f}  additive_gap={gap:+.5f}  "
          f"{'TIGHT' if abs(gap)<1e-3 else 'LOOSE (integrality gap > 0)'}")
    print(f"    rank(X*)={rk} / {n+1}   top eigs {np.round(eig[:5],3)}")
    print(f"    truly-pinned vars (const over ALL ground states): {len(pinned)} / {n}")
    print(f"    SDP |corr|>0.999 candidates: {len(sdp_cand)};  verified-removable: {len(verified)}")
    if pinned:
        kinds = {}
        for i in pinned:
            k = Q.kind[i]; kinds[k] = kinds.get(k, 0) + 1
        print(f"    pinned var kinds: {kinds}")
    return dict(s=s, p=model['p'], n=n, ngs=len(states), sdp=val, gap=gap,
                rank=rk, npin=len(pinned), nverif=len(verified))


if __name__ == "__main__":
    rows = [run(s) for s in (4, 6, 8)]
    print("\n--- modmul integrality-gap vs size ---")
    print(f"{'s':>2} {'p':>4} {'vars':>5} {'#gs':>7} {'SDP_opt':>10} {'gap':>9} {'rankX':>6} {'pinned':>7} {'removable':>9}")
    for r in rows:
        print(f"{r['s']:>2} {r['p']:>4} {r['n']:>5} {r['ngs']:>7} {r['sdp']:>10.4f} "
              f"{r['gap']:>9.4f} {r['rank']:>6} {r['npin']:>7} {r['nverif']:>9}")
