"""U22: EXHAUSTIVE route price over every (slot, leaf pair, source) triple.

Every unordered leaf pair has a unique LCA, so the (slot, a, b) triples are exactly the
32,640 leaf pairs.  Two source choices each (carry a's honest point, or b's) -> 65,280
candidates.  Each is built by u20_sweep.price and scored by checker.py's compiled
equations through M's calibrated forward engine.

usage: python3 u22_pairs.py <shard> <nshard>
"""
import sys, time, pickle, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S
import umodel as U

SH = int(sys.argv[1]); NSH = int(sys.argv[2])
jobs = []
for beta in U.SLOTS:
    la, lb = U.tree[beta]
    for a in sorted(U.LIVELEAF[la]):
        for b in sorted(U.LIVELEAF[lb]):
            jobs.append((beta, a, b))
jobs.sort()
print('total (slot,pair) triples: %d ; shard %d/%d' % (len(jobs), SH, NSH))
res = []
t0 = time.time()
mine = [j for k, j in enumerate(jobs) if k % NSH == SH]
for k, (beta, a, b) in enumerate(mine):
    row = [beta, a, b, None, None]
    for i, src in enumerate((a, b)):
        try:
            n, vv, sd = S.price(beta, a, b, src)
            row[3 + i] = n
        except Exception as e:
            row[3 + i] = -1
    res.append(tuple(row))
    if k % 200 == 0:
        el = time.time() - t0
        print('  %5d/%d  %.0fs  eta %.0fs  best=%s' % (k, len(mine), el,
              el / max(k, 1) * (len(mine) - k),
              min([x for r in res for x in r[3:] if x is not None and x >= 0], default=None)))
        sys.stdout.flush()
pickle.dump(res, open('u_pairs_%d.pkl' % SH, 'wb'))
allv = [x for r in res for x in r[3:] if x is not None and x >= 0]
print('DONE shard %d: %d evals, min %d, %.0fs' % (SH, len(allv), min(allv), time.time() - t0))
