#!/usr/bin/env python3
"""Confirm a planted answer was actually FOUND: every exact split of the planted set S
into (beta of size b, alpha of size |S|-b) must appear as a HIT line whose code decodes
to that beta.   usage: ycheckplant.py report.txt  i,j,k,...  b"""
import sys, itertools
rep, Sarg, b = sys.argv[1], sys.argv[2], int(sys.argv[3])
lo = int(sys.argv[4]) if len(sys.argv) > 4 else 0
hi = int(sys.argv[5]) if len(sys.argv) > 5 else 256
S = sorted(int(v) for v in Sarg.split(','))
hits = set()
lines = 0
for line in open(rep):
    t = line.split()
    if t[0] != 'HIT': continue
    lines += 1
    sz = int(t[1]); code = int(t[2])
    hits.add(tuple((code >> (8 * i)) & 0xFF for i in range(sz)))
want = [tuple(c) for c in itertools.combinations(S, b) if lo <= c[0] < hi]
found = [c for c in want if c in hits]
print('planted S       = %s   (|S|=%d)' % (S, len(S)))
print('scan size b     = %d ;  alpha size = %d' % (b, len(S) - b))
print('HIT lines       = %d ; distinct betas = %d' % (lines, len(hits)))
print('exact splits expected = %d, found = %d' % (len(want), len(found)))
for c in want:
    print('   beta=%-28s %s' % (str(list(c)), 'FOUND' if c in hits else 'MISSING'))
ok = len(found) == len(want)
print('\nPLANT TEST: %s' % ('PASS -- the search finds what it is looking for' if ok else 'FAIL'))
sys.exit(0 if ok else 1)
