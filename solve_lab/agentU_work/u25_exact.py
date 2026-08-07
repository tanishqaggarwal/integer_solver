"""U25: the route price at EVERY slot, exhaustive over leaf pairs and over carried values.

Two measured reductions make this exhaustive rather than sampled:
  (A) u23: the price is INDEPENDENT of the value carried on the route -- honest points of
      other leaves, random 296-bit off-curve points, (3,5) and (0,0) all give the identical
      count.  Only matching one of the two leaves' OWN honest point is cheaper (that leaf
      then routes honestly).  So min over all values = min over {src=a, src=b, mixed axes}.
  (B) u24: price(beta,a,b,src=a) depends only on (beta, b) -- the LYING leaf -- not on the
      honest leaf a.  So a slot's |I|*|J| pairs collapse to |I|+|J| distinct prices.

Hence 32,640 (slot,pair) triples x all values  ->  2,114 evaluations, exactly.
usage: python3 u25_exact.py <shard> <nshard>
"""
import sys, time, pickle, collections, random
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S
import umodel as U

SH = int(sys.argv[1]) if len(sys.argv) > 1 else 0
NSH = int(sys.argv[2]) if len(sys.argv) > 2 else 1

jobs = []          # (beta, lying_leaf, honest_leaf)
for beta in U.SLOTS:
    la, lb = U.tree[beta]
    A = sorted(U.LIVELEAF[la]); B = sorted(U.LIVELEAF[lb])
    for b in B:
        jobs.append((beta, b, A[0]))       # src = honest leaf a  -> b lies
    for a in A:
        jobs.append((beta, a, B[0]))       # src = honest leaf b  -> a lies
jobs.sort()
mine = [j for k, j in enumerate(jobs) if k % NSH == SH]
print('jobs %d ; shard %d/%d -> %d' % (len(jobs), SH, NSH, len(mine)))

res = []
t0 = time.time()
for k, (beta, lying, honest) in enumerate(mine):
    la, lb = U.tree[beta]
    if lying in U.LIVELEAF[la]:
        a, b, src = lying, honest, honest      # honest is on the b side, src=b -> a=lying lies
    else:
        a, b, src = honest, lying, honest
    try:
        n, vv, sd = S.price(beta, a, b, src)
    except Exception as e:
        n = -1
    res.append((beta, lying, honest, n))
    if k % 200 == 0:
        el = time.time() - t0
        print('  %4d/%d %.0fs eta %.0fs best=%s' % (k, len(mine), el,
              el / max(k, 1) * (len(mine) - k),
              min([r[3] for r in res if r[3] >= 0], default=None)))
        sys.stdout.flush()
pickle.dump(res, open('u_exact_%d.pkl' % SH, 'wb'))
print('DONE shard %d: %d evals, min %d, %.0fs' % (SH, len(res),
      min(r[3] for r in res if r[3] >= 0), time.time() - t0))
