#!/usr/bin/env python3
"""Link (i) sweep: for SAME-TREE boolean pairs, turn both on and ask the mod-p affine rigidity engine
whether the derived relations are contradictory.  A pair is 'excluded' if the engine reports
conflicts > 0, i.e. the two pin chains force the same wire to two different residues mod p.

This sweep is a SAMPLE, not an exhaustion: there are C(178,2)=15753 same-tree-A pairs and
C(78,2)=3003 same-tree-B pairs, 18756 in total, at ~15 s each.  The script records exactly how many
pairs it has covered so the fraction can be reported honestly.  Checkpointed and resumable.
"""
import sys, os, json, pickle, time, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cfg_rigid2 import build, p

pins = json.load(open(os.path.join(HERE, 'pins.json')))
sup = pickle.load(open(os.path.join(HERE, 'supp.pkl'), 'rb'))
A = sorted(b for b in sup['7715'] if str(b) in pins)
B = sorted(b for b in sup['34554'] if str(b) in pins)
TOTAL_A = len(A) * (len(A) - 1) // 2
TOTAL_B = len(B) * (len(B) - 1) // 2
OUT = os.path.join(HERE, 'sweep_i.json')
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
COORD = {'A': (12186, 16742), 'B': (14853, 24908)}


def pair(j1, j2, tree):
    partner = 5090 if tree == 'A' else 22106
    R = build([j1, j2, partner])
    rel = R['rel']; C = R['CONST']
    vals = []
    for w in COORD[tree]:
        r, a, b = rel(w)
        vals.append(str((a * C[r] + b) % p) if r in C else None)
    return {'tree': tree, 'conflicts': R['conflicts'], 'forced': vals,
            'excluded': R['conflicts'] > 0}


rnd = random.Random(20260807)
plan = []
# adjacent pairs first (consecutive indices in each tree's list), then uniform random pairs
for L, t in ((A, 'A'), (B, 'B')):
    for i in range(len(L) - 1):
        plan.append((L[i], L[i + 1], t))
for L, t, k in ((A, 'A', 4000), (B, 'B', 2000)):
    seen = set()
    while len(seen) < k:
        i, j = rnd.sample(range(len(L)), 2)
        if i > j: i, j = j, i
        seen.add((i, j))
    for i, j in sorted(seen):
        plan.append((L[i], L[j], t))

t0 = time.time()
for j1, j2, tree in plan:
    key = '%d_%d' % (j1, j2)
    if key in res:
        continue
    try:
        r = pair(j1, j2, tree)
    except Exception as ex:
        r = {'tree': tree, 'error': repr(ex)}
    res[key] = r
    json.dump(res, open(OUT, 'w'))
    nA = sum(1 for x in res.values() if x.get('tree') == 'A')
    nB = sum(1 for x in res.values() if x.get('tree') == 'B')
    exc = sum(1 for x in res.values() if x.get('excluded'))
    notexc = [k for k, x in res.items() if x.get('excluded') is False]
    print('pairs done A %d/%d (%.2f%%)  B %d/%d (%.2f%%)  excluded %d  NOT-excluded %d %s  t=%.0fs'
          % (nA, TOTAL_A, 100.0 * nA / TOTAL_A, nB, TOTAL_B, 100.0 * nB / TOTAL_B,
             exc, len(notexc), notexc[:5], time.time() - t0), flush=True)
