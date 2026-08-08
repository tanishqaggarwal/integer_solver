#!/usr/bin/env python3
"""unit_scale.py -- is the modmul failure a budget problem or a wall?

Give simulated annealing 100x the earlier budget on the smallest sizes that
failed, both with the factors clamped (ancilla-filling only) and free.
"""
import random, time, json
from unit_probe import modmul_qubo
from sa import SA

print(f"{'s':>3} {'vars':>6} {'sweeps':>9} {'clamped':>12} {'free':>12} {'sec':>6}")
res = []
for s in (10, 12):
    p, a, b, c, Q, A, Bw = modmul_qubo(s)
    x, _ = Q.witness({}, {}); assert Q.energy(x) == 0
    sa = SA(Q.Q, Q.n)
    clamp = A + Bw; cs = set(clamp)
    x0 = [0] * Q.n
    for t, v in enumerate(A): x0[v] = (a >> t) & 1
    for t, v in enumerate(Bw): x0[v] = (b >> t) & 1
    for sweeps in (100_000, 1_000_000):
        t0 = time.time()
        hc = bc = 10**18; nh = 0
        for r in range(4):
            rr = random.Random(r)
            st = [x0[i] if i in cs else rr.randrange(2) for i in range(Q.n)]
            e, _ = sa.run(sweeps=sweeps, seed=r, x0=st, clamp=clamp)
            hc = min(hc, e); nh += (e == 0)
        nf = 0
        for r in range(4):
            e, _ = sa.run(sweeps=sweeps, seed=90 + r)
            bc = min(bc, e); nf += (e == 0)
        print(f"{s:3d} {Q.n:6d} {sweeps:9,d} {f'{nh}/4 E>={hc}':>12} "
              f"{f'{nf}/4 E>={bc}':>12} {time.time()-t0:6.0f}", flush=True)
        res.append(dict(s=s, sweeps=sweeps, clamped_hits=nh, clamped_best=hc,
                        free_hits=nf, free_best=bc))
        json.dump(res, open('unit_scale.json', 'w'), indent=1)
