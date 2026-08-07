"""Maximum bipartite matching: check atoms <-> free variables occurring linearly in them.
Matched checks are really GATES for their handle; the UNMATCHED checks are the true constraints."""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P

CHECKS = [a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE = set(u for u in range(L.NVARS) if L.definer.get(u) is None)
print(f"checks={len(CHECKS)}  free inputs={len(FREE)}  gates={len(L.definer)}")

# adjacency: check -> free vars occurring LINEARLY (degree exactly 1 in every monomial)
adj = {}
for a in CHECKS:
    cs = []
    for u in L.avars[a]:
        if u not in FREE:
            continue
        if any(mm.count(u) > 1 for mm in L.polys[a]):
            continue
        cs.append(u)
    adj[a] = cs
deg0 = [a for a in CHECKS if not adj[a]]
print(f"checks with NO free var of their own: {len(deg0)}")

# Hopcroft-Karp style augmenting path matching (simple Kuhn's, adequate here)
matchL = {}   # check -> var
matchR = {}   # var -> check
order = sorted(CHECKS, key=lambda a: len(adj[a]))


def try_aug(a, seen):
    for u in adj[a]:
        if u in seen:
            continue
        seen.add(u)
        if u not in matchR or try_aug(matchR[u], seen):
            matchL[a] = u
            matchR[u] = a
            return True
    return False


sys.setrecursionlimit(100000)
t0 = time.time()
for i, a in enumerate(order):
    if a in matchL:
        continue
    try_aug(a, set())
    if i % 2000 == 0:
        print(f"  {i}/{len(order)} matched={len(matchL)} ({time.time()-t0:.0f}s)", flush=True)
print(f"MAXIMUM MATCHING = {len(matchL)} of {len(CHECKS)} checks ({time.time()-t0:.0f}s)")
unm = [a for a in CHECKS if a not in matchL]
print(f"UNMATCHED (true constraints) = {len(unm)}")
print("  first 40:", unm[:40])
json.dump({'matchL': {str(k): v for k, v in matchL.items()}, 'unmatched': unm},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'match.json'), 'w'))
# how many equations do the unmatched checks live in?
eqs = set()
for a in unm:
    eqs |= set(L.atom2eq.get(a, {}))
print(f"unmatched checks touch {len(eqs)} of {L.NEQ} equations")
