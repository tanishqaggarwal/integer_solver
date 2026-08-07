#!/usr/bin/env python3
"""K18: structure probes on the 256 leaf points and the target.
All arithmetic is exact modular integer arithmetic on the constants read out of the file."""
import sys, os, json, time, collections, itertools
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
import fold as FD
P = FD.P
D = FD.points()
L = [(l['sel'], int(l['X']), int(l['Y'])) for l in D['leaves']]
pts = [(x, y) for _, x, y in L]
sel = [s for s, _, _ in L]
T = (int(D['target']['X']), int(D['target']['Y']))
S = {p: i for i, p in enumerate(pts)}
n = len(pts)
print('n leaves', n)

t0 = time.time()
# 1. negations present?
neg = [(i, S[(x, (-y) % P)]) for i, (x, y) in enumerate(pts) if (x, (-y) % P) in S]
print('leaf pairs P_i = -P_j :', len(neg))

# 2. doubling relations 2P_i = P_j
dbl = [(i, S[FD.add(pts[i], pts[i])]) for i in range(n) if FD.add(pts[i], pts[i]) in S]
print('2*P_i lands on another leaf:', len(dbl), dbl[:12])
# 2b. 2P_i = target?
print('2*P_i == T:', [i for i in range(n) if FD.add(pts[i], pts[i]) == T])

# 3. weight-1 / weight-2 hits on T
print('T is a leaf:', T in S)
w2 = []
allpair = {}
for i in range(n):
    for j in range(i + 1, n):
        s2 = FD.add(pts[i], pts[j])
        allpair[s2] = (i, j)
        if s2 == T: w2.append((i, j))
print('pairs summing to T:', w2, '  (%d pairs, %.1fs)' % (len(allpair), time.time() - t0))

# 4. P_i + P_j = P_k  (weight-3 zero relation, up to signs)
hit3 = [(v, S[k]) for k, v in allpair.items() if k in S]
print('P_i+P_j = P_k :', len(hit3), hit3[:10])

# 5. weight-3 hits on T:  P_i+P_j+P_k = T  <=>  T - P_i in pairsums
t1 = time.time()
hit = []
for i in range(n):
    q = FD.add(T, (pts[i][0], (-pts[i][1]) % P))
    if q in allpair: hit.append((i, allpair[q]))
print('triples summing to T:', len(hit), hit[:5], '%.1fs' % (time.time() - t1))

# 6. weight-4 zero relations: P_i+P_j = P_k+P_l  (a collision between two pair-sums)
#    allpair is a dict so equal sums would have collided; count them properly
cnt = collections.Counter()
for i in range(n):
    for j in range(i + 1, n):
        cnt[FD.add(pts[i], pts[j])] += 1
coll = [k for k, c in cnt.items() if c > 1]
print('pair-sum collisions (weight-4 relations):', len(coll))
json.dump({'neg': neg, 'dbl': dbl, 'w2': w2, 'hit3': [[list(a), b] for a, b in hit3]},
          open(K + '/struct.json', 'w'))
print('total %.1fs' % (time.time() - t0))
