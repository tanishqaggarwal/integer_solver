#!/usr/bin/env python3
"""Reconstruct the doubling ladder over the 256 conditional-pin leaves."""
import json, collections, sys
from model import P, S, inv, to_short, load_points, TARGET
from group import add, mul, neg, fit_B

pts, tree, missing, pins = load_points()
S2 = {i: to_short(pts[i]) for i in sorted(pts)}
B = fit_B(list(S2.values()))
idx = {v: k for k, v in S2.items()}
succ = {i: idx[add(S2[i], S2[i])] for i in S2 if add(S2[i], S2[i]) in idx}
starts = [i for i in S2 if i not in set(succ.values())]
pieces = []
for s in starts:
    ch = [s]; seen = {s}
    while ch[-1] in succ and succ[ch[-1]] not in seen:
        ch.append(succ[ch[-1]]); seen.add(ch[-1])
    pieces.append(ch)
pieces.sort(key=len, reverse=True)
sys.stderr.write('pieces %s\n' % [len(c) for c in pieces])

# splice: end of piece a, doubled twice (one unknown point in between), = start of piece b
startpt = {S2[p[0]]: k for k, p in enumerate(pieces)}
nxt = {}
gapval = {}
for k, pc in enumerate(pieces):
    e = S2[pc[-1]]
    g = add(e, e)              # the missing leaf
    gg = add(g, g)
    if gg in startpt:
        nxt[k] = startpt[gg]; gapval[k] = g
sys.stderr.write('splices %d of %d\n' % (len(nxt), len(pieces)))
head = [k for k in range(len(pieces)) if k not in set(nxt.values())]
sys.stderr.write('heads %s\n' % head)
order = []
if len(head) == 1:
    k = head[0]
    while True:
        order.append(k)
        if k in nxt: k = nxt[k]
        else: break
sys.stderr.write('spliced chain covers %d pieces\n' % len(order))

ladder = []            # list of (exponent, boolean_id_or_None, point)
for j, k in enumerate(order):
    for b in pieces[k]:
        ladder.append((b, S2[b]))
    if k in nxt:
        ladder.append((None, gapval[k]))
sys.stderr.write('ladder length %d  (named %d, unnamed %d)\n'
                 % (len(ladder), sum(1 for b, _ in ladder if b is not None),
                    sum(1 for b, _ in ladder if b is None)))
if __name__ == '__main__':
    json.dump({'B': str(B),
               'ladder': [[b, str(pt[0]), str(pt[1])] for b, pt in ladder]},
              open('ladder.json', 'w'))
    print('wrote ladder.json, length', len(ladder))
