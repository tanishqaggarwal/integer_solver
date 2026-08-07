"""STRUCTURAL SOLVE ORDER.

31,475 gate atoms each define one variable; 7,273 variables are free; 10,792 atoms are checks.
A check can only be *satisfied* by determining something, so pair each check with a free
variable occurring linearly in it (maximum bipartite matching).  Matched checks then act as
extra 'gates'; unmatched checks must hold automatically.  Whether the resulting definition
graph is acyclic decides whether the instance solves by a single topological pass.
"""
import sys, os, collections, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
CHECKS=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
print(f"checks={len(CHECKS)} free vars={len(FREE)} gates={len(L.definer)}")
# bipartite: check -> free vars occurring with degree 1 in every monomial of that atom
adj=collections.defaultdict(list)
for a in CHECKS:
    Pp=L.polys[a]
    cand=set()
    for m in Pp:
        for u in m:
            if u in FREE and m.count(u)==1: cand.add(u)
    for u in cand:
        if all(m.count(u)<=1 for m in Pp):
            adj[a].append(u)
deg=collections.Counter(len(adj[a]) for a in CHECKS)
print("checks by #free-linear candidates:", dict(sorted(deg.items())[:12]), "...")
print("checks with NO free candidate:", deg.get(0,0))
freehit=collections.Counter()
for a in CHECKS:
    for u in adj[a]: freehit[u]+=1
print("free vars reachable from some check:", len(freehit), "of", len(FREE))
# Hopcroft-Karp
t0=time.time()
matchU={}; matchA={}
def bfs():
    dist={}; q=collections.deque()
    for a in CHECKS:
        if a not in matchA: dist[a]=0; q.append(a)
    found=False
    while q:
        a=q.popleft()
        for u in adj[a]:
            b=matchU.get(u)
            if b is None: found=True
            elif b not in dist:
                dist[b]=dist[a]+1; q.append(b)
    return dist, found
def dfs(a, dist, it):
    for u in adj[a]:
        b=matchU.get(u)
        if b is None or (dist.get(b)==dist[a]+1 and dfs(b,dist,it)):
            matchU[u]=a; matchA[a]=u; return True
    dist[a]=-1; return False
sys.setrecursionlimit(100000)
while True:
    dist,found=bfs()
    if not found: break
    prog=False
    for a in CHECKS:
        if a not in matchA and dist.get(a)==0:
            if dfs(a,dist,0): prog=True
    if not prog: break
print(f"maximum matching: {len(matchA)} of {len(CHECKS)} checks ({time.time()-t0:.0f}s)")
print(f"unmatched checks: {len(CHECKS)-len(matchA)}   unused free vars: {len(FREE)-len(matchU)}")
import json
json.dump({'matchA':{str(k):v for k,v in matchA.items()}}, open(os.path.join(HERE,'data','match_struct.json'),'w'))
