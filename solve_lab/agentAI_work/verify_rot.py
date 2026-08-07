#!/usr/bin/env python3
"""Independent check of the rotational-sweep correctness invariant.
Each completed rotation's six shard DONE lines must have candidate counts
summing to exactly C(128,5) = 264566400.  Also checks each shard count
against its own closed form: #{5-subsets of [0,128) with min in [lo,hi)}
 = sum_{m=lo}^{hi-1} C(127-m, 4).
Exit 0 = all completed rotations verified; exit 3 = a violation found."""
import re, sys, os
from math import comb

REP = sys.argv[1] if len(sys.argv) > 1 else \
    '/home/user/integer_solver/solve_lab/agentX_work/rrep_real.txt'
DONEF = os.path.join(os.path.dirname(REP), 'rotdone.txt')

TOTAL = comb(128, 5)
def closed(lo, hi):
    return sum(comb(127 - m, 4) for m in range(lo, min(hi, 128)))

rows = {}
hits = 0
zeroes = 0
pat = re.compile(r'^DONE rot=(\d+) range=\[(\d+),(\d+)\) n=(\d+) zero=(\d+)')
malformed = []
for ln, line in enumerate(open(REP), 1):
    line = line.rstrip('\n')
    if line.startswith('HIT '):
        hits += 1; continue
    m = pat.match(line)
    if not m:
        if line.strip(): malformed.append((ln, line))
        continue
    j, lo, hi, n, z = (int(x) for x in m.groups())
    rows.setdefault(j, {})[(lo, hi)] = n
    zeroes += z

marked = set()
if os.path.exists(DONEF):
    for line in open(DONEF):
        mm = re.match(r'ROT (\d+) DONE', line)
        if mm: marked.add(int(mm.group(1)))

print('C(128,5) = %d' % TOTAL)
print('report file: %s' % REP)
print('rotations with >=1 shard line: %d   marked DONE in rotdone.txt: %d'
      % (len(rows), len(marked)))
print('HIT lines: %d   total zero= field: %d' % (hits, zeroes))
if malformed:
    print('MALFORMED LINES: %r' % malformed[:5])
print()

bad = []
verified = set()
for j in sorted(rows):
    sh = rows[j]
    s = sum(sh.values())
    ok_sum = (s == TOTAL)
    ok_n = (len(sh) == 6)
    percol = []
    for (lo, hi), n in sorted(sh.items()):
        c = closed(lo, hi)
        percol.append(((lo, hi), n, c, n == c))
    ok_each = all(p[3] for p in percol)
    status = 'OK' if (ok_sum and ok_n and ok_each) else 'VIOLATION'
    flag = ' [marked DONE]' if j in marked else ' [not marked done]'
    print('rot %3d: shards=%d sum=%d  delta=%+d  per-shard-closed-form=%s  %s%s'
          % (j, len(sh), s, s - TOTAL,
             'all match' if ok_each else 'MISMATCH', status, flag))
    if not ok_each:
        for p in percol:
            if not p[3]:
                print('     range=%s reported n=%d  closed form=%d  delta=%+d'
                      % (p[0], p[1], p[2], p[1] - p[2]))
    if status == 'OK':
        verified.add(j)
    else:
        bad.append(j)

print()
# rotations marked DONE but not arithmetically verified
unbacked = sorted(marked - verified)
print('rotations marked DONE AND invariant-verified : %d  %s'
      % (len(verified & marked), sorted(verified & marked)))
print('rotations marked DONE but NOT verified       : %d  %s' % (len(unbacked), unbacked))
print('rotations with data but invariant VIOLATED   : %d  %s' % (len(bad), bad))
sys.exit(3 if (bad or unbacked) else 0)
