#!/usr/bin/env python3
"""For every free variable that moves the three EC constraints, measure how many
OTHER constraints it breaks (mod p).  Cheap knobs = candidates for a real solve."""
import os, pickle, random, sys, time
import jengine as E, jman as J, jmodp as MP
from collections import deque

P = MP.P
definer = J.definer
EC = [20407, 20409, 31575]
CONS = MP.CONS

base = [x % P for x in J.BASE]
MP.fwd_modp(base)
r0 = MP.residues(base)
bad0 = set(i for i, x in r0.items() if x)
print("base violated:", sorted(bad0))

# cone leaves of the EC constraints
seen = set(); q = deque()
for i in EC:
    q.extend(E.varsof[i])
leaves = set()
while q:
    x = q.popleft()
    if x in seen: continue
    seen.add(x)
    d = definer.get(x)
    if d is None:
        leaves.add(x); continue
    for w in E.varsof[d]:
        if w != x and w not in seen: q.append(w)

random.seed(2024)
res = []
t0 = time.time()
LV = sorted(leaves)
for n, z in enumerate(LV):
    v = [x % P for x in J.BASE]
    v[z] = random.randrange(P)
    MP.fwd_modp(v)
    r = MP.residues(v)
    moved_ec = [i for i in EC if r[i] != r0[i]]
    if not moved_ec:
        continue
    broke = sorted(i for i in CONS if r[i] and i not in bad0)
    healed = sorted(i for i in bad0 if not r[i])
    res.append((z, len(broke), broke[:8], healed))
    if n % 50 == 0:
        print(f"  {n}/{len(LV)} {time.time()-t0:.0f}s", file=sys.stderr)
res.sort(key=lambda t: t[1])
print(f"{len(res)} movers, sorted by #constraints broken:")
for z, nb, br, hl in res[:60]:
    print(f"  x_{z}: breaks {nb} {br} heals {hl}")
pickle.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'jknobs.pkl'), 'wb'))
