#!/usr/bin/env python3
"""Exact fractional coverage of the c0 m<=8 run, from engine evidence only.

A chunk counts as covered only when ALL EIGHT shard passes for it carry the engine's own
DONE line with the exact expected candidate count for that [lo,hi) range.  Shell markers are
not consulted.  Partial chunks are reported separately and are NOT counted as coverage.
"""
import json, os, glob
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
CH = json.load(open(os.path.join(HERE, 'chunks4.json')))
TOT = comb(256, 4) * 2**4
cov = 0; part = []; hits = 0; zero = 0; secs = 0.0; done_units = 0
for j, (lo, hi) in enumerate(CH):
    n = sum(comb(255 - (s >> 1), 3) * 8 for s in range(lo, hi))
    ok = 0
    for s in range(8):
        f = os.path.join(HERE, 'runs8', 'c0m8.c%d.s%d.txt' % (j, s))
        if not os.path.exists(f): continue
        for ln in open(f):
            p = ln.split()
            if p and p[0] == 'DONE':
                kv = dict(q.split('=', 1) for q in p if '=' in q)
                if kv.get('range') == '[%d,%d)' % (lo, hi) and int(kv['n']) == n:
                    ok += 1; zero += int(kv['zero']); secs += float(p[-1].rstrip('s'))
            if p and p[0] == 'HIT': hits += 1
    done_units += ok
    if ok == 8: cov += n
    elif ok: part.append((j, ok, n))
    print('chunk %2d range=[%3d,%3d) n_each=%-11d shards_exact=%d/8 %s'
          % (j, lo, hi, n, ok, 'COVERED' if ok == 8 else ''))
print('\nCOVERED (all 8 shards): %d of %d candidates = %.3f%% of the m=8 space'
      % (cov, TOT, 100.0 * cov / TOT))
print('partially-run chunks (NOT counted): %s' % (part or 'none'))
print('shard-units complete: %d of %d   engine seconds: %.0f   dx=0 events: %d   HITS: %d'
      % (done_units, 8 * len(CH), secs, zero, hits))
