"""Agent B: maximum bipartite matching gate -> output variable.
Gives the exact number of ASSERTION gates (unmatched) and FREE inputs (unmatched vars)."""
import pickle, collections, sys
sys.setrecursionlimit(100000)

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs = M['facs']
NV = 38748

cand = []
for p in facs:
    quad = set()
    for m in p:
        if len(m) > 1: quad.update(m)
    c = [m[0] for m, co in p.items() if len(m) == 1 and abs(co) == 1 and m[0] not in quad]
    cand.append(c)
print("gates:", len(facs), " total candidate edges:", sum(len(c) for c in cand))
print("gates with 0 candidates:", sum(1 for c in cand if not c))

# Hopcroft-Karp
INF = float('inf')
matchG = [-1]*len(facs)     # gate -> var
matchV = [-1]*NV            # var -> gate
import collections as C
def bfs():
    dist = {}
    q = C.deque()
    for g in range(len(facs)):
        if matchG[g] == -1:
            dist[g] = 0; q.append(g)
        else: dist[g] = INF
    found = False
    while q:
        g = q.popleft()
        for v in cand[g]:
            g2 = matchV[v]
            if g2 == -1: found = True
            elif dist[g2] == INF:
                dist[g2] = dist[g]+1; q.append(g2)
    return dist, found

def dfs(g, dist):
    for v in cand[g]:
        g2 = matchV[v]
        if g2 == -1 or (dist.get(g2, INF) == dist[g]+1 and dfs(g2, dist)):
            matchG[g] = v; matchV[v] = g; return True
    dist[g] = INF
    return False

res = 0
while True:
    dist, found = bfs()
    if not found: break
    prog = 0
    for g in range(len(facs)):
        if matchG[g] == -1 and dfs(g, dist): prog += 1
    if prog == 0: break
    res += prog
    print("  matched %d" % res, flush=True)
print("MAX MATCHING:", res)
print("unmatched gates (ASSERTIONS):", len(facs)-res)
print("unmatched vars (FREE INPUTS):", NV-res)
pickle.dump({'matchG': matchG, 'matchV': matchV}, open(W+'match.pkl','wb'), -1)
unm = [g for g in range(len(facs)) if matchG[g] == -1]
print("assertion gate ids (first 60):", unm[:60])
