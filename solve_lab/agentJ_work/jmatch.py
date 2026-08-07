#!/usr/bin/env python3
"""Maximum bipartite matching atoms -> definable variables (Hopcroft-Karp),
then check acyclicity of the induced circuit."""
import os, pickle, sys, time
from collections import deque, defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
polys, defcands = P['polys'], P['defcands']
NA = len(polys); NV = 38748
adj = [list(d) for d in defcands]

INF = float('inf')
matchA = [-1] * NA
matchV = [-1] * NV

def bfs():
    dist = [INF] * NA
    q = deque()
    for a in range(NA):
        if matchA[a] == -1 and adj[a]:
            dist[a] = 0; q.append(a)
    found = False
    while q:
        a = q.popleft()
        for v in adj[a]:
            b = matchV[v]
            if b == -1:
                found = True
            elif dist[b] == INF:
                dist[b] = dist[a] + 1
                q.append(b)
    return dist, found

def dfs(a, dist):
    for v in adj[a]:
        b = matchV[v]
        if b == -1 or (dist[b] == dist[a] + 1 and dfs(b, dist)):
            matchA[a] = v; matchV[v] = a
            return True
    dist[a] = INF
    return False

sys.setrecursionlimit(100000)
t0 = time.time()
res = 0
while True:
    dist, found = bfs()
    if not found:
        break
    prog = 0
    for a in range(NA):
        if matchA[a] == -1 and adj[a]:
            if dfs(a, dist):
                prog += 1
    if prog == 0:
        break
    res += prog
print(f"matching size {res} in {time.time()-t0:.1f}s")
print("unmatched atoms (constraints):", NA - res)
print("unmatched vars (free):", NV - res)
pickle.dump({'matchA': matchA, 'matchV': matchV},
            open(os.path.join(HERE, 'jmatch.pkl'), 'wb'))

# acyclicity: edges other-vars(i) -> matchA[i]
indeg = Counter()
succ = defaultdict(list)
for i in range(NA):
    v = matchA[i]
    if v == -1:
        continue
    vs = set()
    for k in polys[i]:
        vs.update(k)
    for w in vs:
        if w != v:
            succ[w].append(v)
            indeg[v] += 1
q = deque([v for v in range(NV) if indeg[v] == 0])
seen = 0
while q:
    v = q.popleft(); seen += 1
    for w in succ[v]:
        indeg[w] -= 1
        if indeg[w] == 0:
            q.append(w)
print("topologically ordered vars:", seen, "of", NV, "=> cyclic vars:", NV - seen)
