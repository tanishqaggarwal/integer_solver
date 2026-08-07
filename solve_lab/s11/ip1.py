"""IP #1 -- MINIMUM-COST DEFECT PLACEMENT, solved exactly.

    binary y_c = 1  iff constraint c is left violated
    feasibility :  the surviving constraints admit a perfect matching into their controls
                   (equivalently Hall's condition holds on every subset)
    objective   :  | union of the equations occupied by the violated constraints' atoms |

The constraint set of a channel is small (14), so the IP is solved by exact enumeration over
all 2^14 subsets -- this is a certificate, not a heuristic.
"""
import sys, os, json, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
BIT = set(int(b) for b in LD)
S = json.load(open(os.path.join(HERE, 'data', 'last4_nolock.json')))

# constraint -> (atoms that become nonzero if violated, non-bit controls)
GRAPH = {
    'a688':   ([688], [19750]),
    'a1618':  ([1618], [14853]),
    'a29539': ([29539], [14515]),
    'a26731': ([26731], [16742]),
    'a7881':  ([7881], [2751]),
    'a21050': ([21050], [16441]),
    'a26839': ([26839], [18751]),
    'a40065': ([40065], [28955]),
    'a14445': ([14445], None),
    'a27139': ([27139], None),
    'a34580': ([34580], None),
    'a33796': ([33796], None),
    # violating EITHER mirror residual makes the whole group-1 mirror trio nonzero
    'mirror3719':  ([26719, 26721, 26723], None),
    'mirror25118': ([26719, 26721, 26723], None),
}
for k in ['a14445', 'a27139', 'a34580', 'a33796', 'mirror3719', 'mirror25118']:
    GRAPH[k] = (GRAPH[k][0], [u for u in S[k] if u not in BIT])

NODES = sorted(GRAPH)
COST = {}
for k, (atoms, _) in GRAPH.items():
    e = set()
    for a in atoms:
        e |= set(L.atom2eq.get(a, {}))
    COST[k] = e


def max_matching(nodes):
    ml, mr = {}, {}

    def aug(u, seen):
        for w in GRAPH[u][1]:
            if w in seen:
                continue
            seen.add(w)
            if w not in mr or aug(mr[w], seen):
                ml[u] = w
                mr[w] = u
                return True
        return False
    for u in sorted(nodes, key=lambda x: len(GRAPH[x][1])):
        if u not in ml:
            aug(u, set())
    return len(ml)


print(f"constraints: {len(NODES)}")
print(f"full-set matching = {max_matching(NODES)} of {len(NODES)}  -> deficit "
      f"{len(NODES)-max_matching(NODES)}")
print()
best = None
feasible = []
for r in range(0, len(NODES) + 1):
    for sub in itertools.combinations(NODES, r):
        keep = [u for u in NODES if u not in sub]
        if max_matching(keep) != len(keep):
            continue
        eqs = set()
        for u in sub:
            eqs |= COST[u]
        feasible.append((len(eqs), sub))
        if best is None or len(eqs) < best[0]:
            best = (len(eqs), sub)
    if best is not None and r >= (len(NODES) - max_matching(NODES)) + 1:
        break
feasible.sort()
print("cheapest FEASIBLE defect placements (exact IP optimum first):")
for n, sub in feasible[:10]:
    print(f"  cost {n:3d} equations  <- violate {list(sub)}")
print(f"\nIP OPTIMUM for this channel: {best[0]} failing equations "
      f"-> score {L.NEQ - best[0]}   (violate {list(best[1])})")
