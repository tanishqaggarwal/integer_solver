"""Rigorous deficit computation: maximum bipartite matching between the live CONSTRAINTS of
   channel U=0,V=1 (bits 490,91) and their NON-BIT controls (exhaustive scan, s11/last4b.py).
   Bits are excluded because a1430 = x_490^2 - x_490 etc. force every message bit to {0,1}."""
import sys, os, json, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
BIT = set(int(b) for b in LD)
S = json.load(open(os.path.join(HERE, 'data', 'last4_nolock.json')))

# constraints driven earlier in the pipeline (their control is the variable used to drive them)
GRAPH = {
    'a688   (x37892=G0)':   [19750],
    'a1618  (x13682)':      [14853],
    'a29539 (x1308=x14853)': [14515],
    'a26731 (x16742=x19083)': [16742],
    'a7881  (x2751=x1085)': [2751],
    'a21050 (x16441=x4920)': [16441],
    'a26839 (x18751=x33091)': [18751],
    'a40065 (x28955=x11408)': [28955],
}
for k in ['a14445', 'a27139', 'a34580', 'a33796', 'mirror3719', 'mirror25118']:
    GRAPH[k] = [u for u in S[k] if u not in BIT]

nodes = sorted(GRAPH)
print("constraint -> non-bit controls")
for k in nodes:
    print(f"  {k:24s}: {GRAPH[k]}")

matchL, matchR = {}, {}


def aug(u, seen):
    for w in GRAPH[u]:
        if w in seen:
            continue
        seen.add(w)
        if w not in matchR or aug(matchR[w], seen):
            matchL[u] = w
            matchR[w] = u
            return True
    return False


for u in sorted(nodes, key=lambda x: len(GRAPH[x])):
    if u not in matchL:
        aug(u, set())
print(f"\nconstraints = {len(nodes)}   maximum matching = {len(matchL)}   DEFICIT = {len(nodes)-len(matchL)}")
print("unmatched:", [u for u in nodes if u not in matchL])

# Hall violator: minimal set with |N(S)| < |S|
best = None
for r in range(2, len(nodes) + 1):
    for sub in itertools.combinations(nodes, r):
        nb = set()
        for u in sub:
            nb |= set(GRAPH[u])
        if len(nb) < len(sub):
            best = (sub, nb)
            break
    if best:
        break
if best:
    sub, nb = best
    print(f"\nHALL VIOLATOR: {len(sub)} constraints vs {len(nb)} controls")
    print("  constraints:", list(sub))
    print("  controls   :", sorted(nb))
