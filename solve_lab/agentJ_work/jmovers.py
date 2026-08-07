#!/usr/bin/env python3
"""For each of the key reduced constraints, find ALL free variables that move it
(random probe over GF(p)) in the given boolean branch."""
import sys, random
from collections import deque
import jengine as E, jman as J, jmodp as MP, jsolve2 as S

P = MP.P
definer = J.definer
b1 = int(sys.argv[1]) if len(sys.argv) > 1 else 1
b2 = int(sys.argv[2]) if len(sys.argv) > 2 else 1
CS = [int(a) for a in sys.argv[3:]] or [30271, 8583, 22688, 26603, 34370, 27640, 2694,
                                        20407, 20409, 31575, 731, 31571]
val, bad = S.branch(b1, b2)
print(f"branch ({b1},{b2}) violated: {bad}")
random.seed(31)
for c in CS:
    seen = set(); q = deque(E.varsof[c]); lv = set()
    while q:
        x = q.popleft()
        if x in seen: continue
        seen.add(x)
        d = definer.get(x)
        if d is None:
            lv.add(x); continue
        for w in E.varsof[d]:
            if w != x and w not in seen: q.append(w)
    r0 = MP.atom_modp(c, val)
    mv = []
    for z in sorted(lv):
        v2 = list(val); v2[z] = random.randrange(P); MP.fwd_modp(v2)
        if MP.atom_modp(c, v2) != r0:
            mv.append(z)
    print(f"a{c}: cone {len(seen)} leaves {len(lv)} movers {len(mv)}: {mv[:40]}")
