#!/usr/bin/env python3
"""demo.py -- the whole pipeline on cached synthetic instances: recover k, count runs."""
import sys, time
sys.path.insert(0, '.')
from synth.build_lib import get, load
from synth.solve import solve, _residual_dlog

print("SYNTHETIC END-TO-END RECOVERY  (planted key, perfect-annealer oracle per run)")
print(f"{'bits':>5} {'mu/run':>7} {'runs to hit':>12} {'worst 2^(b-mu)':>15} {'k recovered = planted':>22} {'sec':>6}")
have = sorted(int(b) for b in load())
for b in have:
    inst = get(b)
    mu = min(16, b)
    r = solve(inst, mu=mu, order='planted_first')
    ok = r['found'] and inst.curve.mul(r['k'], inst.G) == inst.T and r['k'] == inst.k
    print(f"{b:5d} {mu:7d} {r['runs']:12d} {('2^%d'%(b-mu)):>15} {'YES' if ok else 'NO':>22} {r['secs']:6.2f}")

print()
print("SOLUTION COUNT: over ALL 2^(b-mu) prefixes, 1 or 2 yield a valid k (k and k+n)")
for b in [bb for bb in have if bb <= 28]:
    inst = get(b); mu = min(12, b); c,G,n,T,k = inst.curve,inst.G,inst.n,inst.T,inst.k
    hits = 0; rec = None
    for khi in range(1 << (b - mu)):
        Tres = c.add(T, c.mul((-(khi<<mu)) % n, G))
        x = _residual_dlog(c, G, n, Tres, mu)
        if x is not None and c.mul(((khi<<mu)+x)%n, G) == T: hits += 1; rec = (khi<<mu)+x
    print(f"  bits={b:3d} mu={mu:2d}: prefixes={1<<(b-mu):6d}  successful={hits}  "
          f"recovered={rec}  planted={k}  match={rec==k}")
